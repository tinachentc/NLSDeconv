import scanpy as sc
import pandas as pd
import torch

from preprocess import Preprocess
from deconv import Deconv
from plot import separate_plt, overall_plt

if __name__ == '__main__':
    # assign scRNA and ST directories (assign celltype_key in scRNA data)
    sc_dir = './example_data/sc.h5ad'
    st_dir = './example_data/st_10000.h5ad'
    celltype_key = 'celltype'

    ad_sc = sc.read_h5ad(sc_dir)
    ad_st = sc.read_h5ad(st_dir)

    # preprocess
    ad_st, ad_sc = Preprocess(ad_st, ad_sc, celltype_key=celltype_key).preprocess()

    # deconvolution
    flist = ['astrocytes', 'Olig', 'endo_mural', 'iNeuron', 'eNeuron', 'microglia']
    # method of Softthreshold Least Square
    res, time_res, head_res = Deconv(ad_sc, ad_st, flist=flist).SLS()
    # method of Nonnegative Least Square (estimate lr when lr=None)
    # for CPU
    # res, time_res, head_res = Deconv(ad_sc, ad_st, flist=flist, normalization=True).NLS(reg=1e-1, lr=1e-2, warm_start=True, num_epochs=1000, device="cpu")
    # for GPU (e.g. cuda:0)
    # res, time_res, head_res = Deconv(ad_sc, ad_st, flist=flist, normalization=True).NLS(reg=1e-1, lr=1e-2, warm_start=True, num_epochs=1000, device="cuda:0")

    # plot
    res = res.cpu().numpy()
    separate_plt(res, head_res, ad_st, show_head_res=head_res, spot_size=400) #separate_plt(res, head_res, ad_st, show_head_res=['iNeuron', 'eNeuron'], spot_size=400)
    overall_plt(res, head_res, ad_st, spot_size=250, margin_size=400)

    # (optional) calculate RMSE
    true_result_path = './example_data/Out_cell_ratio_1x.csv'
    true_result = pd.read_csv(true_result_path, header=0, index_col=0)
    true_result = torch.tensor(true_result.dropna().to_numpy())
    d1, d2 = true_result.shape
    rmse = (true_result - res) ** 2
    rmse = torch.sqrt(rmse.sum() / d1 / d2)
    print(f'RMSE: {rmse.item()}')