# NLSDeconv: Cell-type Deconvolution

## Introduction:
The provided codes implement NLSDeconv, a novel cell-type deconvolution method for spatial transcriptomics (ST) data. 
We have two algorithm options: soft-thresholding least squares (SLS) and non-negative least squares (NLS). SLS is developed as a fast approximation version of NLS, and is recommended for users without GPU access.
To implement the algorithm, user needs to provide an ST dataset and a reference scRNA-seq dataset (with cell type information).

We provide example codes for running SLS on a seqFISH+ dataset in `main_example.py`.

## Requirements:
Developed on `python = 3.11` `PyTorch = 2.0.1`
You can install the rest of the requirements via
`pip install -r requirements.txt`

## Pipeline:
1. Load ST and scRNA-seq datasets separately.
```bash
ad_st = sc.read_h5ad(st_dir)
ad_sc = sc.read_h5ad(sc_dir)
```
2. Preprocess both datasets with one line code, which includes:
 - total count normalization for scRNA-seq dataset
 - removal of cell types with small number of observations for scRNA-seq dataset
 - selection of top genes characterizing each cell type for scRNA-seq dataset
 - matching genes for ST and scRNA-seq datasets
```bash
ad_st, ad_sc = Preprocess(ad_st, ad_sc, celltype_key='celltype').preprocess()
```
3. Apply a deconvolution algorithm, SLS or NLS.
 - SLS
 ```bash
 res, time_res, head_res = Deconv(ad_sc, ad_st).SLS()
 ```
 - NLS
 ```bash
 res, time_res, head_res = Deconv(ad_sc, ad_st, normalization=True).NLS(reg=1e-1, lr=1e-2, warm_start=True, num_epochs=1000)
 ```
4. Visulize the deconvolution result.
 - Separate figures of each cell type
 - An overall pie plot

## Arguments:


## Attribution:

