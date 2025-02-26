# import libraries
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import pickle
import rpy2.robjects as robjects
robjects.r('library(Matrix)')

# preprocess scATAC-seq data using pycisTopic
from pycisTopic.cistopic_class import *
count_matrix=pd.read_csv('pbmccount.tsv', sep='\t')

# creating a cisTopic object(removed blacklist regions)
path_to_blacklist='hg38-blacklist.v2.bed'
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, path_to_blacklist=path_to_blacklist)

# adding metadata to a cisTopic object
cell_data =  pd.read_csv('celltype.tsv', sep='\t')
cistopic_obj.add_cell_data(cell_data)

# save a cisTopic object
work_dir="scenicplus"
if not os.path.exists(os.path.join(work_dir, 'scATAC')):
    os.makedirs(os.path.join(work_dir, 'scATAC'))


pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


# load a cisTopic object
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

# run models
models=run_cgs_models(cistopic_obj,
                    n_topics=[2,4,10,16,32,48],
                    n_cpu=5,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = '/picb/bigdata/project/liangxuan/temopory')

# save models results
if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb'))


# model selection
models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'rb'))
from pycisTopic.lda_models import *
model = evaluate_models(models,
                       select_model=16,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

# clustering and visualization
from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['celltype'])

# topic binarization 
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

# Differentially Accessible Regions (DARs)
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions, split_pattern = '-')

# save region sets
if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))

import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))


# generate a custom cisTarget database

# load pycisTopic region sets
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))

# prepare region sets
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))



# motif collection(motif can be downloaded from: https://resources.aertslab.org/cistarget/databases)
db_fpath = "/picb/bigdata/project/liangxuan"
motif_annot_fpath = "/picb/bigdata/project/liangxuan"
rankings_db = os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')

# create cistarget databases
if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = os.path.join('/picb/bigdata/project/liangxuan/temopory', 'ray_spill'),
    annotation_version = 'v10nr_clust',
    )




# preprocess scRNA-seq data using Scanpy
# import libraries
import pandas as pd
import scanpy as sc
import pandas as pd

# load data
adata = sc.read_h5ad('pbmc.h5ad')

# QC
sc.pp.filter_cells(adata,min_genes=200)
sc.pp.filter_genes(adata, min_cells=0)
import scrublet
sc.external.pp.scrublet(adata)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# visualization
sc.pl.umap(adata,color ='celltype')

# select highly variable features
gene=pd.read_csv("pbmc_gene.csv")
gene1 = gene['x']
gene2= ', '.join(gene1.astype(str))
gene3=pd.Index(gene2.split(', '))
all_genes=adata.var_names
gene_in_known_list = []

for gene in all_genes:
    if gene in gene3:
        gene_in_known_list.append(True)
    else:
        gene_in_known_list.append(False)

adata = adata[:, gene_in_known_list]

# save data
adata.write( 'pbmc1.h5ad', compression='gzip')


# run SCENIC+ 
# load data
import dill
adata = sc.read_h5ad('pbmc1.h5ad')
menr = dill.load(open('scenicplus/motifs/menr.pkl', 'rb'))
cistopic_obj = dill.load(open('scenicplus/scATAC/cistopic_obj.pkl', 'rb'))

# create a scenicplus object
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
    bc_transform_func = lambda x: f'{x}___cisTopic'
) 

scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

# run scenicplus
from scenicplus.wrappers.run_scenicplus import run_scenicplus
from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS
get_search_space(scplus_obj,
                 biomart_host = 'http://www.ensembl.org',
                 species = 'hsapiens',
                 assembly = 'hg38',
                 upstream = [1000, 250000],
                 downstream = [1000, 250000],
                 extend_tss=[1000, 1000],
                 remove_promoters=True)
calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = 15,
                    _temp_dir = '/picb/bigdata/project/liangxuan/temopory',
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)

df_interact=scplus_obj.uns['region_to_gene']

df_interact.to_csv("scenicplus_pbmc_peak_gene.rda", index=False)