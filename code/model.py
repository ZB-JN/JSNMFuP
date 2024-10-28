from tkinter import W
import torch  
import numpy as np
import scanpy as sc
import scvi
from anndata import AnnData
from utils import *
import time 
from scipy.sparse import isspmatrix
# import bct
from bct.algorithms.modularity import *
import random
import umap.umap_ as umap
import matplotlib.pyplot as plt
import os
import pandas as pd


torch.set_default_dtype(torch.float64)





class JSNMFuP:
	'''
	JSNMF model.

	Parameters
	----------
	RNA        : AnnData object that contains the data from RNA, true cell type should be
			      stored in RNA.obs['celltype']
	ATAC       : AnnData object that contains the data from ATAC, or other modality

	max_epochs : Number of max epochs to run, optinal, default is 200 epochs

	'''

	def __init__(self,
				 RNA: AnnData,
				 ATAC: AnnData,
				 R,
				 k: int = 20,
				 max_epochs: int = 200):

		self.RNA = RNA
		self.ATAC = ATAC
		self.R = R
		self.label = RNA.obs['celltype'].values

		self.max_epochs = max_epochs
		self.num_c  =k
		self.embedding = None
		try:
			if isspmatrix(RNA.X):
				self.X1 = torch.from_numpy(RNA.X.toarray()).T.cuda()
				self.X2 = torch.from_numpy(ATAC.X.toarray()).T.cuda()
			else:
				self.X1 = torch.from_numpy(RNA.X).T.cuda()
				self.X2 = torch.from_numpy(ATAC.X).T.cuda()

		except:
			if isspmatrix(RNA.X):
				self.X1 = torch.from_numpy(RNA.X.toarray()).T
				self.X2 = torch.from_numpy(ATAC.X.toarray()).T
				self.R = torch.from_numpy(R.toarray())
			else:
				self.X1 = torch.from_numpy(RNA.X).T
				self.X2 = torch.from_numpy(ATAC.X).T
				self.R = torch.from_numpy(R.toarray())

	def parameter_selection(self):
		'''
		Calculates the initilization values for the model.

		Returns
		----------
		alpha   : init hyperparameter alpha
		beta    : init hyperparameter beta
		Inits   : dict stroing the init values

		'''
		Inits = {}
		num_c = self.num_c
		W1, H1 = nndsvd(self.X1, num_c)
		print('nnsvd W1 done')
		W2, H2 = nndsvd(self.X2, num_c)
		print('nnsvd W2 done')

		Inits['W1'] = W1
		Inits['W2'] = W2
		Inits['H1'] = H1
		Inits['H2'] = H2

		D1 = dist2(self.X1.T,self.X1.T)
		print('D1 done')
		S1 = affinityMatrix(D1,20) # cell by cell
		print('S1 done')
		D2 = dist2(self.X2.T,self.X2.T)
		print('D2 done')
		S2 = affinityMatrix(D2,20) #cell by cell
		print('S2 done')

		Inits['A1'] = S1
		Inits['D1'] = torch.diag(S1.sum(0))
		Inits['L1'] = Inits['D1'] - Inits['A1']
		Inits['A2'] = S2
		Inits['D2'] = torch.diag(S2.sum(0))
		Inits['L2'] = Inits['D2'] - Inits['A2']

		print('SNF starts')
		W = SNF(S1,S2,K = 20)
		print('SNF done')
		Inits['S'] = W

		H1tH1 = H1.T @ H1
		H2tH2 = H2.T @ H2

		err1 = torch.square(torch.norm(self.X1 - W1 @ H1))
		err2 = torch.square(torch.norm(self.X2 - W2 @ H2))
		err3 = 0.5 * (torch.square(torch.norm(Inits['S'] - H1tH1)) + torch.square(torch.norm(Inits['S'] - H2tH2)))
		err4 = 0.5 * (H1 @ Inits['L1'] @ H1.T + H2 @ Inits['L2'] @ H2.T).trace()
		err5 = (W1.T @ self.R @ W2).trace()
		alpha = (err1 + err2) / err3
		gamma = torch.abs((err1 + err2) / err4)
		gamma2 = torch.abs((err1 + err2) / err5)
		alpha = alpha / 10
		gamma = gamma / 10
		gamma2 = gamma2 / 10

		return alpha, gamma, gamma2,  Inits

	def gcnmf(self):
		'''
		Main function to run the model.

		Returns
		----------
        W1,W2,H1,H2    : resulting values of init  after running the model
		S			   : resulting complete graph

		use_epochs     : used epochs to converge

		objs		   : value of objective during each iteration

		'''
		Maxiter = self.max_epochs
		X1 = self.X1
		X2 = self.X2
		R = self.R
		alpha = self.alpha
		gamma = self.gamma
		gamma2 = self.gamma2
		W1 = self.Inits['W1']
		W2 = self.Inits['W2']
		H1 = self.Inits['H1'].T
		H2 = self.Inits['H2'].T
		A1 = self.Inits['A1']
		D1 = self.Inits['D1']
		S = self.Inits['S']
		S = S / S.sum(dim=0, keepdim=True)
		A2 = self.Inits['A2']
		D2 = self.Inits['D2']
		n = S.size(0)

		try:
			theta1 = torch.tensor(1 / 2).cuda()
			theta2 = torch.tensor(1 / 2).cuda()
			yita1 = torch.tensor(1.).cuda()
			yita2 = torch.tensor(0.5).cuda()
			Lamda = torch.tenosor(0.5).cuda()
			obj_old = torch.tensor(1.).cuda()
			objs = torch.zeros((Maxiter, 1)).cuda()
		except:
			theta1 = torch.tensor(1 / 2)
			theta2 = torch.tensor(1 / 2)
			yita1 = torch.tensor(1.)
			yita2 = torch.tensor(0.5)
			Lamda = torch.tensor(0.5)
			obj_old = torch.tensor(1.)
			objs = torch.zeros((Maxiter, 1))

		for i in range(Maxiter):
			# update W1 W2 using bpp algorithm
			# W1, _ = nnlsm_blockpivot(H1, X1.T, False, W1.T)
			# W2, _ = nnlsm_blockpivot(H2, X2.T, False, W2.T)
			# W1 = W1.T
			# W2 = W2.T

			# updating w1 w2 via multiplication rule
			H1tH1 = H1.T @ H1
			H2tH2 = H2.T @ H2
			tmp1 = W1 @ H1tH1 + 0.5 * gamma2 * R @ W2
			tmp1[tmp1 < 1e-10] = 1e-10
			W1 = W1 * (X1 @ H1) / tmp1
			tmp2 = W2 @ H2tH2 + 0.5 * gamma2 * R.T @ W1
			tmp2[tmp2 < 1e-10] = 1e-10
			W2 = W2 * (X2 @ H2) / tmp2

			# update H1, H2
			W1tW1 = W1.T @ W1
			W2tW2 = W2.T @ W2
			H1tH1 = H1.T @ H1
			H2tH2 = H2.T @ H2
			tmp_deno_1 = gamma * theta1.square() * D1 @ H1 + (alpha + yita1) * H1 @ H1tH1 + H1 @ W1tW1
			tmp_deno_1[tmp_deno_1 < 1e-10] = 1e-10
			tmp_nume_1 = alpha * S.T @ H1 + X1.T @ W1 + gamma * theta1.square() * A1 @ H1 + yita1 * H1
			H1 = H1 * (tmp_nume_1 / tmp_deno_1)

			tmp_deno_2 = gamma * theta2.square() * D2 @ H2 + (alpha + yita2) * H2 @ H2tH2 + H2 @ W2tW2
			tmp_deno_2[tmp_deno_2 < 1e-10] = 1e-10
			tmp_nume_2 = alpha * S.T @ H2 + X2.T @ W2 + gamma * theta2.square() * A2 @ H2 + yita2 * H2
			H2 = H2 * (tmp_nume_2 / tmp_deno_2)

			# update S
			H1tH1 = H1 @ H1.T
			H2tH2 = H2 @ H2.T
			Q = 1 / 2 * (H1tH1 + H2tH2)
			tmp = alpha*S + Lamda * S.sum(dim = 0, keepdim =True)
			tmp[tmp<1e-10] = 1e-10
			try:
				tmp_ones = torch.ones((n,n)).cuda()
			except:
				tmp_ones = torch.ones((n,n))
			S = S * ((alpha * Q + Lamda * tmp_ones) / tmp)


			# update theta
			tmp1 = 1 / torch.trace(H1.T @ (D1 - A1) @ H1)
			tmp2 = 1 / torch.trace(H2.T @ (D2 - A2) @ H2)
			tmp = 1 / torch.trace(H1.T @ (D1 - A1) @ H1) + 1 / torch.trace(H2.T @ (D2 - A2) @ H2)
			theta1 = tmp1 / tmp
			theta2 = tmp2 / tmp

			# if stop_rule == 2:
			obj = compute_obj(X1, X2, R, S, W1, W2, H1, H2, D1, A1, D2, A2, alpha, gamma, gamma2, theta1,
							  theta2)
			objs[i, 0] = obj
			error = torch.abs(obj_old - obj) / obj_old
			if (error < 1e-5 and i > 0) or i == Maxiter - 1:
			# if i == Maxiter - 1:
				print('number of epoch:', i + 1)
				print('obj:', obj)
				print('converged!')
				break

			obj_old = obj

			print('number of epoch:', i + 1)
			print('obj:', obj)

		S = (S + S.T) / 2
		use_epochs = i + 1

		return W1, W2, H1, H2, S, use_epochs, objs

	def run(self):
		'''
		Run the JSNMF model. Init time and main function time are recorded.

		Returns
		----------
		result  :  dict storing some information during the model running

		'''
		start = time.time()
		alpha, gamma, gamma2, Inits = self.parameter_selection()
		end = time.time()
		self.init_t = end - start
		self.alpha = alpha
		self.gamma = gamma
		self.gamma2 = gamma2
		self.Inits = Inits

		print('Init done')

		start = time.time()
		W1, W2, H1, H2, S, used_epoch, objs = self.gcnmf()
		end = time.time()
		self.run_t = end - start
		self.used_epoch = used_epoch

		result = dict(W1=W1, W2=W2,
					  H1=H1, H2=H2,
					  S=S, used_epoch=used_epoch,
					  objs=objs,
					  init_t=self.init_t,
					  run_t=self.run_t)

		self.result = result

	# return result

	def cluster(self, K=50, step=0.01, start=2.7, upper=4, seed=3):

		'''
		Use louvain to cluster the cells based on the complete graph S, note that
		the function tries to find  a partition that has the same number of clusters
		as the true labels, the resolution parameter of louvain is found using binary search

		Parameters
		----------
		K      : (0, N) int, parameter of Wtrim
			     Number of neighbors to retain
		step   :  the step of binary search to find the partition, default 0.01
		start  :  start searching point of binary search
		upper  :  the upper bound of the reolution paramter to be searched
		seed   : seed parameter for louvain algorithm

		Returns
		-------
		res_clu : resulting cluster labels for each cell

		Note that sometimes exact number of clusters as true labels may not be found, paramters
		need to be adjusted then, like step, seed and upper

		'''
		S = self.result['S']
		A = Wtrim(S, K=50)
		A = A.numpy()

		# num_c = self.num_c
		num_c= len(np.unique(self.label))

		tmp_gamma = start
		# use louvain to cluster based on A
		clusters, q_stat = community_louvain(A, gamma=tmp_gamma, seed=seed)
		tmp_c = len(np.unique(clusters))
		tmp_clu = clusters

		res_clu = None

		# use binary search to find the corret gamma parameter
		while True:
			if tmp_c == num_c:
				res_clu = tmp_clu
				break

			if tmp_c < num_c:
				tmp_gamma = tmp_gamma + step
			else:
				tmp_gamma = tmp_gamma - step

			if tmp_gamma < 0 or tmp_gamma > upper:
				break
			tmp_gamma = round(tmp_gamma, 2)
			# print(tmp_res)
			clusters, q_stat = community_louvain(A, gamma=tmp_gamma, seed=seed)
			tmp_c = len(np.unique(clusters))
			tmp_clu = clusters

		return res_clu

	def v_embedding(self, label, tag=False, **kwarg):

		'''
		Visualize based on the complete graph S using Umap

		Parameters
		----------
		label     : array, true or clustered (louvain result) labels for each cell
		tag		  : if recalculte umap embedding
		**kwarg   : kwarg for the umap

		Returns
		-------
		res_clu : resulting cluster labels for each cell

		'''

		# transfer S to distance matrix first
		S = self.result['S']
		data = 1 - S
		data = data - np.diag(data.diagonal())
		reducer = umap.UMAP(**kwarg)

		# avoid recompute umap embedding
		if self.embedding is None:
			# min_dist = 0.68, n_neighbors=12
			embedding = reducer.fit_transform(data)
			self.embedding = embedding


	def enrich_analysis(self, genes=None, peaks=None, topk_gene=200, topk_peak=1000,
						folder='h3k4'):
		'''
		Parameters
		----------
		genes     : gene names, if not none, save the factors for gene enrichment analysis
		peaks	  : peaks names, if not none, save the factors for region enrichment analysis
		topk_gene : number of top genes for each factor, default 200
		topk_peak : number of top peak for each factor, default 1000
		folder    : the name of the folder to save the results

		Returns
		-------
		saved files containing factors for

		'''

		if genes is not None:
			path = folder + '/W1'
			if not os.path.exists(path):
				os.makedirs(path)
			W1 = self.result['W1']
			nfactor = W1.shape[1]
			sorted, indices = torch.sort(W1, dim=0, descending=True)
			print('Start writing factors for gene enrichment analysis')
			for i in range(nfactor):
				tmp_dict = {'gene symbol': genes[indices[:topk_gene, i]]}
				tmp_pd = pd.DataFrame(tmp_dict)
				tmp_pd.to_csv(folder + '/W1/factor_' + str(i + 1) + '.csv', index=False)
				print('Successfully write factor_' + str(i + 1) + ' for gene enrichment analysis')
			print('done')
			print()

		if peaks is not None:
			path = folder + '/W2'
			if not os.path.exists(path):
				os.makedirs(path)
			W2 = self.result['W2']
			nfactor = W2.shape[1]
			sorted, indices = torch.sort(W2, dim=0, descending=True)
			print('Start writing factors for region enrichment analysis')
			for i in range(nfactor):
				tmp_dict = {'peak loc': peaks[indices[:topk_peak, i]]}
				tmp_pd = pd.DataFrame(tmp_dict)
				tmp_pd.to_csv(folder + '/W2/factor_' + str(i + 1) + '.csv')
				print('Successfully write factor_' + str(i + 1) + ' for region enrichment analysis')

			print('done')






            
class JSNMFuP_3mod:

	def __init__(self,
				 X1: AnnData,
				 X2: AnnData,
				 X3: AnnData,
				 label,R12,R13,R23,
				 num_c: int =20,
				 max_epochs: int = 200):

		if isspmatrix(X1.X):
			self.X1 = torch.from_numpy(X1.X.toarray()).T
			self.X2 = torch.from_numpy(X2.X.toarray()).T
			self.X3 = torch.from_numpy(X3.X.toarray()).T
			self.R12 = torch.from_numpy(R12.toarray())
			self.R13 = torch.from_numpy(R13.toarray())
			self.R23 = torch.from_numpy(R23.toarray())
		else:
			self.X1 = torch.from_numpy(X1.X).T
			self.X2 = torch.from_numpy(X2.X).T
			self.X3 = torch.from_numpy(X3.X).T
			self.R12 = torch.from_numpy(R12.toarray())
			self.R13 = torch.from_numpy(R13.toarray())
			self.R23 = torch.from_numpy(R23.toarray())

		self.label = label

		self.max_epochs = max_epochs
		# self.num_c = num_c = len(np.unique(self.label))
		self.num_c = num_c
		self.embedding = None

	def parameter_selection(self):
		'''
		Calculates the initilization values for the model.

		Returns
		----------
		alpha   : init hyperparameter alpha
		beta    : init hyperparameter beta
		Inits   : dict stroing the init values

		'''
		Inits = {}
		num_c = self.num_c
		W1, H1 = nndsvd(self.X1, num_c)
		print('nnsvd W1 done')
		W2, H2 = nndsvd(self.X2, num_c)
		print('nnsvd W2 done')
		W3, H3 = nndsvd(self.X3, num_c)
		print('nnsvd W3 done')

		Inits['W1'] = W1
		Inits['W2'] = W2
		Inits['W3'] = W3
		Inits['H1'] = H1
		Inits['H2'] = H2
		Inits['H3'] = H3



		D1 = dist2(self.X1.T,self.X1.T)
		print('D1 done')
		S1 = affinityMatrix(D1,20) # cell by cell
		print('S1 done')

		D2 = dist2(self.X2.T,self.X2.T)
		print('D2 done')
		S2 = affinityMatrix(D2,20) #cell by cell
		print('S2 done')

		D3 = dist2(self.X3.T,self.X3.T)
		print('D3 done')
		S3 = affinityMatrix(D3,20) #cell by cell
		print('S3 done')



		Inits['A1'] = S1
		Inits['D1'] = torch.diag(S1.sum(0))
		Inits['L1'] = Inits['D1'] - Inits['A1']

		Inits['A2'] = S2
		Inits['D2'] = torch.diag(S2.sum(0))
		Inits['L2'] = Inits['D2'] - Inits['A2']

		Inits['A3'] = S3
		Inits['D3'] = torch.diag(S3.sum(0))
		Inits['L3'] = Inits['D3'] - Inits['A3']



		print('SNF starts')
		W = SNF(S1,S2,S3,K = 20)
		# W = SNF(Inits['A1'], Inits['A2'], Inits['A3'],K=20)
		print('SNF done')
		Inits['S'] = W

		H1tH1 = H1.T @ H1
		H2tH2 = H2.T @ H2
		H3tH3 = H3.T @ H3

		err1 = torch.square(torch.norm(self.X1 - W1 @ H1))
		err2 = torch.square(torch.norm(self.X2 - W2 @ H2))
		err21 = torch.square(torch.norm(self.X3 - W3 @ H3))
		err3 = 1 / 3 * (torch.square(torch.norm(Inits['S'] - H1tH1)) + torch.square(
			torch.norm(Inits['S'] - H2tH2)) + torch.square(torch.norm(Inits['S'] - H3tH3)))
		err4 = 1 / 3 * (H1 @ Inits['L1'] @ H1.T + H2 @ Inits['L2'] @ H2.T + H3 @ Inits['L3'] @ H3.T).trace()
		err5= (W1.T@self.R12@W2).trace()
		err6 = (W1.T @ self.R13 @ W3).trace()
		err7 = (W2.T @ self.R23 @ W3).trace()
		alpha = (err1 + err2 + err21) / err3
		gamma = torch.abs((err1 + err2 + err21) / err4)
		# beta1=torch.abs((err1 + err2 + err21) / err5)
		# beta2 = torch.abs((err1 + err2 + err21) / err6)
		# beta3 = torch.abs((err1 + err2 + err21) / err7)
		beta1 = torch.abs(3*(err1 + err2 + err21) / (err5+err6+err7))
		beta3 =beta2 =beta1
		alpha = alpha / 10
		gamma = gamma / 10
		beta1 =beta1/10
		beta2=beta2/10
		beta3=beta3/10

		return alpha, gamma, beta1, beta2, beta3, Inits


	def gcnmf(self):
		'''
		Main function to run the model.

		Returns
		----------
        W1,W2,H1,H2    : resulting values of init  after running the model
		S			   : resulting complete graph

		use_epochs     : used epochs to converge

		objs		   : value of objective during each iteration

		'''
		Maxiter = self.max_epochs
		X1 = self.X1
		X2 = self.X2
		X3 = self.X3
		R12 = self.R12
		R13 = self.R13
		R23 = self.R23
		alpha = self.alpha
		gamma = self.gamma
		beta1 = self.beta1
		beta2 = self.beta2
		beta3 = self.beta3
		W1 = self.Inits['W1']
		W2 = self.Inits['W2']
		W3 = self.Inits['W3']
		H1 = self.Inits['H1'].T
		H2 = self.Inits['H2'].T
		H3 = self.Inits['H3'].T
		A1 = self.Inits['A1']
		D1 = self.Inits['D1']
		S = self.Inits['S']
		S = S / S.sum(dim=0, keepdim=True)
		A2 = self.Inits['A2']
		D2 = self.Inits['D2']
		A3 = self.Inits['A3']
		D3 = self.Inits['D3']

		n = S.size(0)

		try:
			theta1 = torch.tensor(1 / 2).cuda()
			theta2 = torch.tensor(1 / 2).cuda()
			theta3 = torch.tensor(1 / 3).cuda()
			yita1 = torch.tensor(1.).cuda()
			yita2 = torch.tensor(0.5).cuda()
			yita3 = torch.tensor(0.5).cuda()
			Lamda = torch.tenosor(0.5).cuda()
			obj_old = torch.tensor(1.).cuda()
			objs = torch.zeros((Maxiter, 1)).cuda()
		except:
			theta1 = torch.tensor(1 / 3)
			theta2 = torch.tensor(1 / 3)
			theta3 = torch.tensor(1 / 3)
			yita1 = torch.tensor(1.)
			yita2 = torch.tensor(0.5)
			yita3 = torch.tensor(0.5)
			Lamda = torch.tensor(0.5)
			obj_old = torch.tensor(1.)
			objs = torch.zeros((Maxiter, 1))

		for i in range(Maxiter):
			# update W1 W2 using bpp algorithm
			# W1, _ = nnlsm_blockpivot(H1, X1.T, False, W1.T)
			# W2, _ = nnlsm_blockpivot(H2, X2.T, False, W2.T)
			# W3, _ = nnlsm_blockpivot(H3, X3.T, False, W3.T)
			# W1 = W1.T
			# W2 = W2.T
			# W3 = W3.T

			# updating w1 w2 via multiplication rule
			H1tH1 = H1.T @ H1
			H2tH2 = H2.T @ H2
			H3tH3 = H3.T @ H3
			tmp1 = W1 @ H1tH1+0.5*beta1*R12@W2+0.5*beta2*R13@W3
			tmp1[tmp1 < 1e-10] = 1e-10
			W1 = W1 * (X1 @ H1) / tmp1
			tmp2 = W2 @ H2tH2+0.5*beta1*R12.T@W1+0.5*beta3*R23@W3
			tmp2[tmp2 < 1e-10] = 1e-10
			W2 = W2 * (X2 @ H2) / tmp2
			tmp3 = W3 @ H3tH3+0.5*beta2*R13.T@W1+0.5*beta3*R23.T@W2
			tmp3[tmp3 < 1e-10] = 1e-10
			W3 = W3 * (X3 @ H3) / tmp3

			# update H1, H2
			W1tW1 = W1.T @ W1
			W2tW2 = W2.T @ W2
			W3tW3 = W3.T @ W3
			H1tH1 = H1.T @ H1
			H2tH2 = H2.T @ H2
			H3tH3 = H3.T @ H3

			tmp_deno_1 = gamma * theta1.square() * D1 @ H1 + (alpha + yita1) * H1 @ H1tH1 + H1 @ W1tW1
			tmp_deno_1[tmp_deno_1 < 1e-10] = 1e-10
			tmp_nume_1 = alpha * S.T @ H1 + X1.T @ W1 + gamma * theta1.square() * A1 @ H1 + yita1 * H1
			H1 = H1 * (tmp_nume_1 / tmp_deno_1)

			tmp_deno_2 = gamma * theta2.square() * D2 @ H2 + (alpha + yita2) * H2 @ H2tH2 + H2 @ W2tW2
			tmp_deno_2[tmp_deno_2 < 1e-10] = 1e-10
			tmp_nume_2 = alpha * S.T @ H2 + X2.T @ W2 + gamma * theta2.square() * A2 @ H2 + yita2 * H2
			H2 = H2 * (tmp_nume_2 / tmp_deno_2)

			tmp_deno_3 = gamma * theta3.square() * D3 @ H3 + (alpha + yita3) * H3 @ H3tH3 + H3 @ W3tW3
			tmp_deno_3[tmp_deno_3 < 1e-10] = 1e-10
			tmp_nume_3 = alpha * S.T @ H3 + X3.T @ W3 + gamma * theta3.square() * A3 @ H3 + yita3 * H3
			H3 = H3 * (tmp_nume_3 / tmp_deno_3)

			# update S
			H1tH1 = H1 @ H1.T
			H2tH2 = H2 @ H2.T
			H3tH3 = H3 @ H3.T
			Q = 1 / 3 * (H1tH1 + H2tH2 + H3tH3)
			tmp = alpha * S + Lamda * S.sum(dim=0, keepdim=True)
			tmp[tmp < 1e-10] = 1e-10
			try:
				tmp_ones = torch.ones((n, n)).cuda()
			except:
				tmp_ones = torch.ones((n, n))
			S = S * ((alpha * Q + Lamda * tmp_ones) / tmp)

			# update theta
			tmp1 = 1 / torch.trace(H1.T @ (D1 - A1) @ H1)
			tmp2 = 1 / torch.trace(H2.T @ (D2 - A2) @ H2)
			tmp3 = 1 / torch.trace(H3.T @ (D3 - A3) @ H3)
			tmp = tmp1 + tmp2 + tmp3
			theta1 = tmp1 / tmp
			theta2 = tmp2 / tmp
			theta3 = tmp3 / tmp

			# if stop_rule == 2:
			obj = compute_obj_3_R(X1, X2, X3, R12,R13,R23,S, W1, W2, W3, H1, H2, H3, D1, A1, D2, A2, D3, A3, alpha, gamma, beta1,beta2,beta3,theta1,
								theta2, theta3)
			objs[i, 0] = obj
			error = torch.abs(obj_old - obj) / obj_old
			if (error < 1e-5 and i > 0) or i == Maxiter - 1:
				print('number of epoch:', i + 1)
				print('obj:', obj)
				print('converged!')
				break

			obj_old = obj

			print('number of epoch:', i + 1)
			print('obj:', obj)

		S = (S + S.T) / 2
		use_epochs = i + 1

		return W1, W2, W3, H1, H2, H3, S, use_epochs, objs

	def run(self):
		'''
		Run the JSNMF model. Init time and main function time are recorded.

		Returns
		----------
		result  :  dict storing some information during the model running

		'''
		start = time.time()
		alpha, gamma, beta1,beta2,beta3,Inits = self.parameter_selection()
		# alpha, gamma, beta1, beta2, Inits = self.parameter_selection()
		end = time.time()
		self.init_t = end - start
		self.alpha = alpha
		self.gamma = gamma
		self.beta1 = beta1
		self.beta2 = beta2
		self.beta3 = beta3
		self.Inits = Inits

		print('Init done')

		start = time.time()
		W1, W2, W3, H1, H2, H3, S, use_epochs, objs = self.gcnmf()
		end = time.time()
		self.run_t = end - start
		self.used_epoch = use_epochs

		result = dict(W1=W1, W2=W2, W3=W3,
					  H1=H1, H2=H2,
					  H3=H3,
					  S=S, used_epoch=use_epochs,
					  objs=objs,
					  init_t=self.init_t,
					  run_t=self.run_t)

		self.result = result

	# return result

	def cluster(self, K=50, step=0.01, start=2.7, upper=4, seed=3):

		'''
		Use louvain to cluster the cells based on the complete graph S, note that
		the function tries to find  a partition that has the same number of clusters
		as the true labels, the resolution parameter of louvain is found using binary search

		Parameters
		----------
		K      : (0, N) int, parameter of Wtrim
			     Number of neighbors to retain
		step   :  the step of binary search to find the partition, default 0.01
		start  :  start searching point of binary search
		upper  :  the upper bound of the reolution paramter to be searched
		seed   : seed parameter for louvain algorithm

		Returns
		-------
		res_clu : resulting cluster labels for each cell

		Note that sometimes exact number of clusters as true labels may not be found, paramters
		need to be adjusted then, like step, seed and upper

		'''
		S = self.result['S']
		A = Wtrim(S, K=50)
		A = A.numpy()

		# num_c = self.num_c
		num_c = len(np.unique(self.label))
		tmp_gamma = start
		# use louvain to cluster based on A
		clusters, q_stat = community_louvain(A, gamma=tmp_gamma, seed=seed)
		tmp_c = len(np.unique(clusters))
		tmp_clu = clusters


		return tmp_clu

	def v_embedding(self, label, tag=False, **kwarg):

		'''
		Visualize based on the complete graph S using Umap

		Parameters
		----------
		label     : array, true or clustered (louvain result) labels for each cell
		tag		  : if recalculte umap embedding
		**kwarg   : kwarg for the umap

		Returns
		-------
		res_clu : resulting cluster labels for each cell

		'''

		# transfer S to distance matrix first
		S = self.result['S']
		data = 1 - S
		data = data - np.diag(data.diagonal())
		reducer = umap.UMAP(**kwarg)

		# avoid recompute umap embedding
		if self.embedding is None:
			# min_dist = 0.68, n_neighbors=12
			embedding = reducer.fit_transform(data)
			self.embedding = embedding

		# recaculate embedding if needed
		if tag is True:
			embedding = reducer.fit_transform(data)
			self.embedding = embedding


	def enrich_analysis(self, genes=None, peaks=None, topk_gene=200, topk_peak=200,
						folder='sctrio'):
		'''
		Parameters
		----------
		genes     : gene names, if not none, save the factors for gene enrichment analysis
		peaks	  : peaks names, if not none, save the factors for region enrichment analysis
		topk_gene : number of top genes for each factor, default 200
		topk_peak : number of top peak for each factor, default 1000
		folder    : the name of the folder to save the results

		Returns
		-------
		saved files containing factors for

		'''

		if genes is not None:
			path = folder + '/W1'
			if not os.path.exists(path):
				os.makedirs(path)
			W1 = self.result['W1']
			nfactor = W1.shape[1]
			sorted, indices = torch.sort(W1, dim=0, descending=True)
			print('Start writing factors for gene enrichment analysis')
			for i in range(nfactor):
				tmp_dict = {'gene symbol': genes[indices[:topk_gene, i]]}
				tmp_pd = pd.DataFrame(tmp_dict)
				tmp_pd.to_csv(folder + '/W1/factor_' + str(i + 1) + '.csv', index=False)
				print('Successfully write factor_' + str(i + 1) + ' for gene enrichment analysis')
			print('done')
			print()

		if peaks is not None:
			path = folder + '/W2'
			if not os.path.exists(path):
				os.makedirs(path)
			W2 = self.result['W2']
			nfactor = W2.shape[1]
			sorted, indices = torch.sort(W2, dim=0, descending=True)
			print('Start writing factors for region enrichment analysis')
			for i in range(nfactor):
				tmp_dict = {'peak loc': peaks[indices[:topk_peak, i]]}
				tmp_pd = pd.DataFrame(tmp_dict)
				tmp_pd.to_csv(folder + '/W2/factor_' + str(i + 1) + '.csv')
				print('Successfully write factor_' + str(i + 1) + ' for region enrichment analysis')

			print('done')



