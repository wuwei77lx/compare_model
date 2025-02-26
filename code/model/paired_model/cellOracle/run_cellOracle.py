# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
co.__version__

# load scATAC peak data and peak connection data made with Cicero

# load scATAC peak data
peaks = pd.read_csv("cellOracle_pbmc_all_peaks.csv", index_col=0)
peaks = peaks.x.values
peaks

# load peak connection data
cicero_connections = pd.read_csv("cellOracle_pbmc_cicero_connections.csv", index_col=0)
cicero_connections.head()

# annotate transcription start sites(TSS)
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38") 
tss_annotated

# Integrate TSS info and cicero connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, 
                                               cicero_connections=cicero_connections)
integrated.to_csv("cellOracle_pbmc_peak_gene.csv")