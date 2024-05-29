import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

from anndata import AnnData
from scipy import sparse

class baseline(nn.Module):
    def __init__(self, init_mat, input_dim, output_dim):
        super(baseline, self).__init__()
        self.weight = nn.Parameter(torch.Tensor(output_dim, input_dim))
        if init_mat is None:
            nn.init.xavier_uniform_(self.weight)
        else:
            assert (output_dim, input_dim) == init_mat.shape, "initialization dimensions don't match the required ones"
            with torch.no_grad():
                self.weight.data = init_mat

    def forward(self, value):
        output = torch.matmul(self.weight, value)
        return output

def loss_function(yhat, y):
    return torch.norm(yhat-y, p='fro')/2

class Deconv:
    def __init__(self, ad_sc, ad_st, celltype_key='celltype', flist=None, normalization=False):
        assert isinstance(ad_sc, AnnData), "ad_sc is not an instance of AnnData"
        assert isinstance(ad_st, AnnData), "ad_st is not an instance of AnnData"
        assert celltype_key in ad_sc.obs.columns, "celltype_key not found in ad_sc.obs"
        self.sc_count = ad_sc.X
        self.st_count = ad_st.X
        self.input_dim, self.embed_dim = self.sc_count.shape
        self.output_dim, embed_dim2 = self.st_count.shape
        assert self.embed_dim == embed_dim2, "column dimensions of datasets don't match"

        self.celltype_mapping = ad_sc.obs[celltype_key]
        self.num_celltypes = len(self.celltype_mapping.unique())

        if flist is None:
            self.flist = list(self.celltype_mapping.unique())
        else:
            #flist = [x.replace('.', '_') for x in flist]##this is needed for survey paper data
            self.flist = flist
        if sparse.isspmatrix(self.sc_count):
            self.sc_count = self.sc_count.todense()
        if sparse.isspmatrix(self.st_count):
            self.st_count = self.st_count.todense()

        self.m = None
        if normalization:
            self.m = max(self.sc_count.max(), self.st_count.max())
            self.sc_count = self.sc_count / self.m
            self.st_count = self.st_count / self.m

    def SLS(self):
        start_LS = time.time()

        beta, _, _, _ = np.linalg.lstsq(self.sc_count.T, self.st_count.T, rcond=None)
        M_LS = torch.tensor(beta.T, dtype=torch.float64)

        result_LS = torch.zeros(self.output_dim, self.num_celltypes)
        for i in range(self.num_celltypes):
            c = self.flist[i]
            columns = [j for j, cell_type in enumerate(self.celltype_mapping) if cell_type == c]
            result_LS[:, i] = M_LS[:, columns].sum(dim=1)
        result_LS.clamp_(min=0)
        res_LS = result_LS / result_LS.sum(dim=1, keepdim=True)
        res_LS[torch.isnan(res_LS)] = 1 / self.num_celltypes

        end_LS = time.time()
        time_LS = end_LS - start_LS
        print(f'Time for SLS: {time_LS}')

        return res_LS, time_LS, self.flist

    def NLS(self, lr, reg=0., warm_start=True, num_epochs=5000, device="cpu"):
        if device != "cpu" and torch.cuda.is_available():
            print(f'Using device: {device}')
            device = torch.device(device)
        else:
            print(f'Using device: cpu')
            device = torch.device("cpu")

        start_GD = time.time()

        # learning rate preparation
        if lr is None:
            singular_values = np.linalg.svd(self.sc_count, compute_uv=False)
            lr = 1/singular_values[0] ** 2
        print(f'Using lr: {lr:.4f}')

        # init preparation
        if warm_start:
            beta, _, _, _ = np.linalg.lstsq(self.sc_count.T, self.st_count.T, rcond=None)
            init_mat = torch.tensor(beta.T, dtype=torch.float64).float()
            init_mat = F.relu(init_mat)
            if self.m is not None:
                init_mat /= self.m
        else:
            init_mat = None

        # data preparation
        input_dat = torch.tensor(self.sc_count, dtype=torch.float64)
        input_dat = input_dat.to(device).float()
        output_dat = torch.tensor(self.st_count, dtype=torch.float64)
        output_dat = output_dat.to(device).float()

        # model setting
        model0 = baseline(init_mat=init_mat, input_dim=self.input_dim, output_dim=self.output_dim)
        model0.to(device)
        optimizer = optim.Adam(model0.parameters(), lr=lr)#optim.SGD(model0.parameters(), lr=lr)

        # training
        for epoch in range(num_epochs):
            optimizer.zero_grad()
            output = model0(input_dat)
            loss = loss_function(output, output_dat)+reg*torch.norm(model0.weight, p='fro')/2
            loss.backward()
            optimizer.step()
            print(f'Epoch {epoch + 1}, Loss: {loss.item():.4f}')

            with torch.no_grad():
                model0.weight.data.clamp_(min=0)

        # result
        model0.eval()
        with torch.no_grad():
            res0 = model0.weight.data

            res_GD = torch.zeros(self.output_dim, self.num_celltypes)

            for i in range(self.num_celltypes):
                c = self.flist[i]
                columns = [j for j, cell_type in enumerate(self.celltype_mapping) if cell_type == c]
                res_GD[:, i] = res0[:, columns].sum(dim=1)

            res_GD = res_GD / res_GD.sum(dim=1, keepdim=True)
            res_GD[torch.isnan(res_GD)] = 1 / self.num_celltypes

        end_GD = time.time()
        time_GD = end_GD - start_GD
        print(f'Time for NLS: {time_GD}')

        return res_GD, time_GD, self.flist
