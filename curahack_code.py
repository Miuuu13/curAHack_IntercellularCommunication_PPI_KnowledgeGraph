#%%
# 
#imports (pre-installed libraries)
import pandas as pd
#import numpy as np
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




####### loading and filtering .csv from cellchat completed #######

#%%

""" Visualize the PPI interaction in form of a network """




# Build a directed graph from ligand→receptor, labeled by source→target
G = nx.DiGraph()

for _, row in merged_up_df.iterrows():
    ligand = row['ligand']
    receptor = row['receptor']
    source = row['source']
    target = row['target']
    label = f"{source} → {target}"
    G.add_edge(ligand, receptor, label=label)

# Plot it
plt.figure(figsize=(14, 10))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=800, font_size=10)
plt.title("Ligand–Receptor Network")
plt.axis('off')
plt.show()


#%%

""" Improved graph - UP """

# step 1, filtering already done 


# -----------------------------
# Step 2: Build NetworkX Graph
# -----------------------------
G = nx.DiGraph()

for _, row in merged_up_df.iterrows():
    ligand = row['ligand']
    receptor = row['receptor']
    prob = row['prob']
    
    # Add edge with 'prob' as weight
    G.add_edge(ligand, receptor, weight=prob)

# -----------------------------
# Step 3: Visualization Settings
# -----------------------------

# Node sizes: based on degree
node_sizes = [G.degree(n) * 100 for n in G.nodes()]

# Edge widths: scaled by 'prob'
edge_widths = [d['weight'] * 20 for _, _, d in G.edges(data=True)]

# Layout
pos = nx.spring_layout(G, seed=42)

# Plot
plt.figure(figsize=(18, 14))
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='red', alpha=0.9)
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color='black', arrows=True)
nx.draw_networkx_labels(G, pos, font_size=9)

plt.title("Ligand–Receptor Network (Filtered, Weighted by Prob)", fontsize=16)
plt.axis('off')
plt.tight_layout()

# Create folder if it doesn't exist
os.makedirs("figures", exist_ok=True)

# Get current timestamp (e.g., 2025-03-31_14-22-05)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# Construct file paths
svg_path = f"figures/ligand_receptor_network_up_{timestamp}.svg"
png_path = f"figures/ligand_receptor_network_up_{timestamp}.png"

# Save the plots
plt.savefig(svg_path, format='svg')
plt.savefig(png_path, format='png', dpi=300)

# Show the plot
plt.show()

#save figure
plt.savefig("ligand_receptor_network.svg", format='svg')
plt.savefig("ligand_receptor_network.png", format='png')
# Show the plot
plt.show()


#%%



""" Improved graph - DOWN """


# step 1, filtering already done 


# -----------------------------
# Step 2: Build NetworkX Graph
# -----------------------------
G = nx.DiGraph()

for _, row in merged_down_df.iterrows():
    ligand = row['ligand']
    receptor = row['receptor']
    prob = row['prob']
    
    # Add edge with 'prob' as weight
    G.add_edge(ligand, receptor, weight=prob)

# -----------------------------
# Step 3: Visualization Settings
# -----------------------------

# Node sizes: based on degree
node_sizes = [G.degree(n) * 100 for n in G.nodes()]

# Edge widths: scaled by 'prob'
edge_widths = [d['weight'] * 20 for _, _, d in G.edges(data=True)]

# Layout
pos = nx.spring_layout(G, seed=42)

# Plot
plt.figure(figsize=(18, 14))
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='skyblue', alpha=0.9)
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color='black', arrows=True)
nx.draw_networkx_labels(G, pos, font_size=9)

plt.title("Ligand–Receptor Network (Filtered, Weighted by Prob)", fontsize=16)
plt.axis('off')
plt.tight_layout()

# Create folder if it doesn't exist
os.makedirs("figures", exist_ok=True)

# Get current timestamp (e.g., 2025-03-31_14-22-05)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# Construct file paths
svg_path = f"figures/ligand_receptor_network_down_{timestamp}.svg"
png_path = f"figures/ligand_receptor_network_down_{timestamp}.png"

# Save the plots
plt.savefig(svg_path, format='svg')
plt.savefig(png_path, format='png', dpi=300)

# Show the plot
plt.show()

#save figure
plt.savefig("ligand_receptor_network.svg", format='svg')
plt.savefig("ligand_receptor_network.png", format='png')
# Show the plot
plt.show()




#%%

""" Check out hihly up/down regulated using  logFC"""

# set thresholds

# Thresholds for logFC (tweak as needed)
UP_LOGFC_THRESHOLD = 0.5
DOWN_LOGFC_THRESHOLD = -0.5
# then filter
# --- Highly upregulated interactions ---

highly_upregulated_df = merged_up_df[
    (merged_up_df['ligand.logFC'] >= UP_LOGFC_THRESHOLD) |
    (merged_up_df['receptor.logFC'] >= UP_LOGFC_THRESHOLD)
].reset_index(drop=True)

# --- Highly downregulated interactions ---
highly_downregulated_df = merged_down_df[
    (merged_down_df['ligand.logFC'] <= DOWN_LOGFC_THRESHOLD) |
    (merged_down_df['receptor.logFC'] <= DOWN_LOGFC_THRESHOLD)
].reset_index(drop=True)


#%%

""" graph for HIGLY UP"""

# step 1, filtering already done 


# -----------------------------
# Step 2: Build NetworkX Graph
# -----------------------------
G = nx.DiGraph()

for _, row in highly_upregulated_df.iterrows():
    ligand = row['ligand']
    receptor = row['receptor']
    prob = row['prob']
    
    # Add edge with 'prob' as weight
    G.add_edge(ligand, receptor, weight=prob)

# -----------------------------
# Step 3: Visualization Settings
# -----------------------------

# Node sizes: based on degree
node_sizes = [G.degree(n) * 100 for n in G.nodes()]

# Edge widths: scaled by 'prob'
edge_widths = [d['weight'] * 20 for _, _, d in G.edges(data=True)]

# Layout
pos = nx.spring_layout(G, seed=42)

# Plot
plt.figure(figsize=(18, 14))
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='red', alpha=0.9)
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color='black', arrows=True)
nx.draw_networkx_labels(G, pos, font_size=9)

plt.title("Ligand–Receptor Network (Filtered, Weighted by Prob) Higly up (logFC > 0.5)", fontsize=16)
plt.axis('off')
plt.tight_layout()

# Create folder if it doesn't exist
os.makedirs("figures", exist_ok=True)

# Get current timestamp (e.g., 2025-03-31_14-22-05)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# Construct file paths
svg_path = f"figures/ligand_receptor_network_down_{timestamp}.svg"
png_path = f"figures/ligand_receptor_network_down_{timestamp}.png"

# Save the plots
plt.savefig(svg_path, format='svg')
plt.savefig(png_path, format='png', dpi=300)

# Show the plot
plt.show()

#save figure
plt.savefig("ligand_receptor_network_HIGHLY_UP.svg", format='svg')
plt.savefig("ligand_receptor_network_HIGHLY_UP.png", format='png')
# Show the plot
plt.show()

#%%

""" graph for HIGLY DOWN """

# step 1, filtering already done 


# -----------------------------
# Step 2: Build NetworkX Graph
# -----------------------------
G = nx.DiGraph()

for _, row in highly_downregulated_df.iterrows():
    ligand = row['ligand']
    receptor = row['receptor']
    prob = row['prob']
    
    # Add edge with 'prob' as weight
    G.add_edge(ligand, receptor, weight=prob)

# -----------------------------
# Step 3: Visualization Settings
# -----------------------------

# Node sizes: based on degree
node_sizes = [G.degree(n) * 100 for n in G.nodes()]

# Edge widths: scaled by 'prob'
edge_widths = [d['weight'] * 20 for _, _, d in G.edges(data=True)]

# Layout
pos = nx.spring_layout(G, seed=42)

# Plot
plt.figure(figsize=(18, 14))
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='skyblue', alpha=0.9)
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color='black', arrows=True)
nx.draw_networkx_labels(G, pos, font_size=9)

plt.title("Ligand–Receptor Network (Filtered, Weighted by Prob) Higly down (low logFC - 0.5)", fontsize=16)
plt.axis('off')
plt.tight_layout()

# Create folder if it doesn't exist
os.makedirs("figures", exist_ok=True)

# Get current timestamp (e.g., 2025-03-31_14-22-05)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# Construct file paths
svg_path = f"figures/ligand_receptor_network_down_{timestamp}.svg"
png_path = f"figures/ligand_receptor_network_down_{timestamp}.png"

# Save the plots
plt.savefig(svg_path, format='svg')
plt.savefig(png_path, format='png', dpi=300)

# Show the plot
plt.show()

#save figure
plt.savefig("ligand_receptor_network_HIGHLY_DOWN.svg", format='svg')
plt.savefig("ligand_receptor_network_HIGHLY_DOWN.png", format='png')
# Show the plot
plt.show()


#%%
################ KNOWLEDGE GRAPH #############################




#%%

""" The already existing knowledge graph """

prime_kg_graph = pd.read_csv("./data/prime_kg.csv")

#%%
prime_kg_graph











#################################################################

# %%
#<ipython-input-8-6617004efd71>:5: DtypeWarning: Columns (3,8) have mixed types. Specify dtype option on import or set low_memory=False.
#  prime_kg_graph = pd.read_csv("./data/prime_kg.csv")


# Define the path to your CSV file
file_path = './data/prime_kg.csv'

# Load the CSV file and specify column data types
prime_kg_df = pd.read_csv(
    file_path,
    dtype={
        'relation': str, 
        'display_relation': str,
        'x_index': float,
        'x_id': str, #int,?
        'x_type': str,
        'x_name': str,
        'x_source': str,
        'y_index': int,
        'y_id': str, #int,?
        'y_type': str,
        'y_name': str,
        'y_source': str,
    }
)

#relation, - protein_protein, 
#display_relation, - ppi,
#x_index, - 0,9796
#x_id, -gene/protein,
#x_type, -PHYHIP,
#x_name, - NCBI
#x_source, 8889
#y_index, 56992
#y_id, -gene/protein
#y_type, -KIF15
#y_name,
#y_source -NCBI

#%%
prime_kg_df


#%%














# %%
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load your data (assuming already loaded as prime_kg_df)
df_subset = prime_kg_df.head(50)  # First 50 rows

# Create the graph
G = nx.Graph()

# Add edges (nodes will be added automatically)
for _, row in df_subset.iterrows():
    G.add_edge(
        row['x_name'], 
        row['y_name'], 
        relation=row['relation'], 
        display_relation=row['display_relation']
    )

# Draw the graph
plt.figure(figsize=(14, 14))
pos = nx.spring_layout(G, seed=42)  # Reproducible layout

# Nodes and edges
nx.draw_networkx_nodes(G, pos, node_size=600, node_color='lightgreen')
nx.draw_networkx_edges(G, pos, alpha=0.4, width=1.5)
nx.draw_networkx_labels(G, pos, font_size=8)

# Optional: draw edge labels
# edge_labels = nx.get_edge_attributes(G, 'display_relation')
# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

plt.title("Protein-Protein Interaction Graph (First 50 rows)", fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.show()

# %%


# %%
