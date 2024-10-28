# JSNMFup_py
JSNMFuP: A unsupervised method for the integrative analysis of single-cell multi-omics data based on non-negative matrix factorization.

This is the Python implementation of the JSNMFuP algorithm. Note that this implementation supports GPU acceleration.

## 1. Installation
You can use the following command to install JSNMFuP:
```
pip install JSNMFuP-py
```

## 2. Usage
The main class "JSNMFuP" needs to be initialized with three and more objects for the "RNA", "ATAC" data and the "the feature association matrix between 'RNA' and 'ATAC' data. The preprocessed data are stored in 'RNA.X', 'ATAC.X' and 'R'. And the real cell labels should be stored in 'RNA.obs['celltype']'. Note that the data preprocessing process is done via python. The maximum number of epochs for a run (i.e., the 'max_epochs' parameter) is set to 200 by default. therefore, it is very simple to initialize the JSNMFuP model with the following code:

```
from JSNMFuP.model import JSNMFuP
test_model = JSNMFuP(rna,atac,R)
```
After initializing, run the model is also quite easy: 
```
test_model.run()
```
The result is saved in 'test_model.result', which is a dict, and the major output of JSNMFuP, the complete graph S, can be get with easy access:
```
S = test_model.result['S']
```
'JSNMFuP' class also has other methods, you can use the 'help' or '?' command for more details explanations of the methods.

