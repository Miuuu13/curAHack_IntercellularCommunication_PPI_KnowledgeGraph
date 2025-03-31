#%%
# 
#imports (pre-installed libraries)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
from datetime import datetime

import networkx as nx
import matplotlib.pyplot as plt

#from pywikigraph import WikiGraph


#%%

# the Cellchat output is seperated in 4 .csv files

down_lig_df = pd.read_csv("./data/net.downLigand_Cellchat.csv")

down_rec_df = pd.read_csv("./data/net.downReceptor_Cellchat.csv")

up_lig_df = pd.read_csv("./data/net.upLigand_Cellchat.csv")

up_rec_df = pd.read_csv("./data/net.upReceptor_Cellchat.csv")


#%%
down_lig_df
#%%
down_rec_df
# %%
up_lig_df
#%%
up_rec_df

#%% 

""" Filter the dfs """

def filter_cellchat_df(df):
    """
    Filters a CellChat dataframe:
    - Keeps only rows where pval <= 0.05
    - Drops rows with missing ligand and receptor percentage values
    """
    return df[
        (df['pval'] <= 0.05) &
        df['ligand.pct.1'].notna() &
        df['ligand.pct.2'].notna() &
        df['receptor.pct.1'].notna() &
        df['receptor.pct.2'].notna()
    ].reset_index(drop=True)

# Apply the filter
down_lig_df = filter_cellchat_df(down_lig_df)
down_rec_df = filter_cellchat_df(down_rec_df)
up_lig_df   = filter_cellchat_df(up_lig_df)
up_rec_df   = filter_cellchat_df(up_rec_df)


#%%
down_lig_df
#%%
down_rec_df
# %%
up_lig_df
#%%
up_rec_df


#%%

""" Merge the two up files """

# up_lig_df and up_rec_df

# Merge (concatenate) the two
merged_up_df = pd.concat([up_lig_df, up_rec_df], ignore_index=True)

# Optional: Drop duplicates if needed
#up_df = up_df.drop_duplicates()

# Preview the result
print(merged_up_df.shape)

print(merged_up_df.head())

#%%

""" Merge the two down files """

# up_lig_df and up_rec_df

# Merge (concatenate) the two
merged_down_df = pd.concat([down_lig_df, down_rec_df], ignore_index=True)

# Optional: Drop duplicates if needed
#up_df = up_df.drop_duplicates()

# Preview the result
print(merged_down_df.shape)

print(merged_down_df.head())
#%%


#%%
# Check size before filtering
print("Merged up:", merged_up_df.shape)
print("Merged down:", merged_down_df.shape)

# Step-by-step inspection
print("Non-NaN ligand.pct.1 in up_df:", merged_up_df['ligand.pct.1'].notna().sum())
print("Non-NaN ligand.pct.2 in up_df:", merged_up_df['ligand.pct.2'].notna().sum())

print("Rows with both pct.1 and pct.2 present:", 
      merged_up_df.dropna(subset=['ligand.pct.1', 'ligand.pct.2']).shape[0])

print("Rows with pval <= 0.05 (no NaN check yet):", 
      merged_up_df[merged_up_df['pval'] <= 0.05].shape[0])


#%%
