import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os

def separate_plt(res, head_res, ad_st, spot_size=400, ncols=3, show_head_res=None, file_name=None):
    # construct plot data
    coords = [[x, y] for x, y in zip(ad_st.obs["x"].values, ad_st.obs["y"].values)]
    adata = sc.AnnData(X=res, obs=coords)
    adata.obsm['spatial'] = np.array(coords)
    adata.var_names = head_res

    # plot
    if show_head_res is None:
        show_head_res = head_res
    sc.pl.spatial(adata, color=show_head_res, cmap="inferno", spot_size=spot_size, vmax=1, vmin=0, ncols=ncols, save=file_name, return_fig=True)

def overall_plt(res, head_res, ad_st, spot_size=250, margin_size=400, file_name=None):
    # construct plot data
    coords = [[x, y] for x, y in zip(ad_st.obs["x"].values, ad_st.obs["y"].values)]
    adata = sc.AnnData(X=res, obs=coords)
    adata.obsm['spatial'] = np.array(coords)
    adata.var_names = head_res

    # plot
    color_map = plt.cm.viridis(np.linspace(0, 1, len(adata.var_names)))

    fig, ax = plt.subplots()
    ax.scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1], s=0)

    # draw a pie chart for each observation
    for idx in range(adata.n_obs):
        ratios = adata.X[idx, :]
        x, y = adata.obsm['spatial'][idx]
        y = max(ad_st.obs["y"].values) + spot_size - y
        #size = spot_size
        #wedges, texts = ax.pie(ratios, colors=plt.cm.viridis(ratios), center=(x, y), radius=np.sqrt(size / np.pi))
        wedges, texts = ax.pie(ratios, colors=plt.cm.viridis(ratios), center=(x, y), radius=spot_size)

    legend_labels = [f'{ct}' for ct in adata.var_names]
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[i]) for i in range(len(legend_labels))]
    ax.legend(handles, legend_labels, loc="lower left", fontsize='small', title_fontsize='medium', prop={'size': 8})

    ax.set_xlim([min(adata.obsm['spatial'][:, 0]) - 10, max(adata.obsm['spatial'][:, 0]) + 10])
    ax.set_ylim([min(adata.obsm['spatial'][:, 1]) - 10 - margin_size, max(adata.obsm['spatial'][:, 1]) + 10])
    #plt.tight_layout()

    if file_name is not None:
        dire = './figures'
        if not os.path.exists(dire):
            os.makedirs(dire)
        plt.savefig(os.path.join(dire, file_name))
    plt.show()