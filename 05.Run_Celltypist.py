import scanpy as sc
import celltypist
from celltypist import models

## set
sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

input_data="./Step1.seurat.data_for_python.h5ad"
plot_dir='./Python.vis/'
data_dir='./'
sc.settings.figdir = plot_dir

input_data="./Step1.seurat.data_for_python.h5ad"
models.models_path 
models.models_description()
model = models.Model.load(model = 'Immune_All_AddPIP.pkl')
model

predictions = celltypist.annotate(adata,
    model = 'COVID19_Immune_Landscape.pkl',
    majority_voting = True)
predictions.predicted_labels.to_csv("{}2_COVID19_Immune_Landscape.csv".format(data_dir))

predictions = celltypist.annotate(adata,
    model = 'Immune_All_Low.pkl',
    majority_voting = True)
predictions.predicted_labels.to_csv("{}1_Immune_All_Low.pkl.csv".format(data_dir))

predictions = celltypist.annotate(adata,
    model = 'Immune_All_AddPIP.pkl',
    majority_voting = True)
predictions.predicted_labels.to_csv("{}3_Immune_All_AddPIP.csv".format(data_dir))

## Our PCCAT annotation model
#predictions = celltypist.annotate(adata,
#    model = 'PCCAT_L2_train_100_model.pkl',
#    majority_voting = True)
#predictions.predicted_labels.to_csv("{}PCCAT_L2_label.csv".format(data_dir))
