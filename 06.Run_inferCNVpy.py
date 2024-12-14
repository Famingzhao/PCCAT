import argparse
import os
from os.path import expanduser

parse=argparse.ArgumentParser(description='inferCNVpy')
parse.add_argument('--sc_input',help='the input file of matrix',type=str,required = True)
parse.add_argument('--sc_format',help = 'the format of input file',choices = ['diopy_h5','h5ad','loom'],required = True)
parse.add_argument('--outdir',help='the output data dir',type=str)
parse.add_argument('--sc_name',help='the name of sc_data',type=str,required = True)
parse.add_argument('--sc_gtf',help='the Cellrange gtf files',type=str,required = True)
parse.add_argument('--reference_key',help='key name for cell types',type=str,required = True)

argv = parse.parse_args()
sc_input = argv.sc_input
sc_format = argv.sc_format
outdir = argv.outdir
sc_name = argv.sc_name
sc_gtf = argv.sc_gtf
reference_key = argv.reference_key

sc_input = os.path.join(expanduser(sc_input))
sc_gtf = os.path.join(expanduser(sc_gtf))
outdir = os.path.join(expanduser(outdir)+"/")
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print(f"Create '{outdir}'")
else:
	print("--- There is this folder! ---")    

############### Run
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
#import diopy
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,figsize=(5, 5))
sc.settings.figdir = outdir

### 1. input data
if sc_format == 'h5':
	adata=sc.read_10x_h5(sc_input,genome=None, gex_only=True)

elif sc_format == 'h5ad':
	adata=sc.read(sc_input)
    
elif sc_format == 'diopy_h5':
	adata=diopy.input.read_h5(sc_input)
 

cnv.io.genomic_position_from_gtf(adata=adata, gtf_file=sc_gtf,
                                 gtf_gene_id='gene_name', inplace=True)

### 2.Run inferCNVpy
## 2.1 We provide all immune cell types as "normal cells".
#np.unique(adata.obs['celltype'])

cnv.tl.infercnv(adata=adata,
    reference_key=reference_key,
    #reference_cat=['B','CD4T', 'CD8T','Erythroid cells','Mast cells', 'Mye', 'NK','Treg', 'Unknown'],
    window_size=100
)

## 2.2 Vis
plt.figure() # creat a new figure
myfig = plt.gcf() # Get the current figure. If no current figure exists, a new one is created using figure().
cnv.pl.chromosome_heatmap(adata, groupby=reference_key,  show=False, dendrogram=True)
plt.savefig("{}Step1.CNV_score_celltype.pdf".format(outdir), dpi=300)

## 2.3 降维聚类
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)

# Vis
plt.figure() # creat a new figure
myfig = plt.gcf() # Get the current figure. If no current figure exists, a new one is created using figure().
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden",show=False, dendrogram=True)
plt.savefig("{}Step2.cnv.pl.chromosome_heatmap_by_cnv_leiden.pdf".format(outdir), dpi=300)

## 2.4 UMAP plot of CNV profiles
cnv.tl.cnv_score(adata)

# Vis
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax3,
    show=False,
)
cnv.pl.umap(adata, color=reference_key, ax=ax1, show=False)
cnv.pl.umap(adata, color="cnv_score", ax=ax4, show=False)
plt.savefig("{}Step3.Vis_cnv_leiden.pdf".format(outdir), dpi=300)

### 3.save results
adata.obs[['cnv_leiden','cnv_score']].to_csv('{}adata_inferCNV.data.csv'.format(outdir))
adata.write(outdir+sc_name+".h5ad")
