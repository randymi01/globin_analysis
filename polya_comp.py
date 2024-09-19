import pandas as pd
import numpy as np


polya_samples = [
    "HS_299",
    "HS_644",
    "HS_672",
    "HS_678",
    "HS_694",
    "HS_696",
    "HS_698",
    "HS_701",
    "HS_702",
    "HS_703"
]

dfs = []


for code in polya_samples:
    path = f"/share/studies/Dermatology_Data/Data_Analysis/HS_Bulk_RNAseq_PBMC/04.Alignment/{code}_1/02_featureCounts2/read_counts_s0_PE.txt.summary"
    df = pd.read_csv(path, delimiter = '\t', names = ["Stat", "Value"])
    df["Sample"] = code+"_pbmc"
    dfs.append(df.copy())
        
merged_df = pd.concat(dfs, ignore_index=True)