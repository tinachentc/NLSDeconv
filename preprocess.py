import scanpy as sc
import numpy as np
from anndata import AnnData
from pandas import DataFrame

class Preprocess:
    def __init__(self, stobj, obj, celltype_key='celltype', cellcount_min=2, gene_top=200, cellcount_norm=True):
        assert isinstance(obj, AnnData), "st is not an instance of AnnData"
        self.stobj = stobj
        assert isinstance(obj, AnnData), "sc is not an instance of AnnData"
        self.obj = obj
        assert celltype_key in obj.obs.columns, "celltype_key not found in sc.obs"
        self.celltype_key = celltype_key
        assert isinstance(cellcount_min, int) and cellcount_min >= 0, "cellcount_min must be a nonnegative integer"
        self.cellcount_min = cellcount_min
        assert isinstance(gene_top, int) and gene_top >= 0, "gene_top must be a nonnegative integer"
        self.gene_top = gene_top
        assert isinstance(cellcount_norm, bool), "cellcount_norm must be a logical variable"
        self.cellcount_norm = cellcount_norm

    def preprocess(self):
        # normalize
        if self.cellcount_norm:
            sc.pp.normalize_total(self.obj)

        # cell cut
        celltype_counts = self.obj.obs[self.celltype_key].value_counts()
        celltype_drop = celltype_counts.index[celltype_counts < self.cellcount_min]
        self.obj = self.obj[~self.obj.obs[self.celltype_key].isin(celltype_drop),].copy()

        # gene cut
        sc.tl.rank_genes_groups(self.obj, groupby=self.celltype_key, use_raw=not self.cellcount_norm)
        markers_df = DataFrame(self.obj.uns["rank_genes_groups"]["names"]).iloc[0:self.gene_top, :]
        genes_sc = np.unique(markers_df.melt().value.values)
        #self.obj = self.obj[:,genes_sc].copy()
        genes = list(set(genes_sc).intersection(set(self.stobj.var_names.values)))
        print(f'Number of genes: {len(genes)}')

        # select common genes
        self.stobj = self.stobj[:, genes].copy()
        self.obj = self.obj[:, genes].copy()

        return self.stobj, self.obj