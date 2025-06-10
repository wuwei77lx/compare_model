# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object

print(co.__version__)  # Print celloracle version to confirm import

# Load scATAC peak data
peaks = pd.read_csv("cellOracle_PBMC_all_peaks.csv", index_col=0)
peaks = peaks.x.values  # Extract peak strings from 'x' column (make sure this column exists)
peaks

# Load peak connection data from Cicero
cicero_connections = pd.read_csv("cellOracle_PBMC_cicero_connections.csv", index_col=0)
cicero_connections.head()

# Annotate transcription start sites (TSS) for peaks
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
tss_annotated.head()

# Integrate TSS annotation with Cicero peak connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, 
                                              cicero_connections=cicero_connections)

# Save integrated peak-gene connections to CSV
integrated.to_csv("CellOracle_PBMC_peak_gene.csv")
