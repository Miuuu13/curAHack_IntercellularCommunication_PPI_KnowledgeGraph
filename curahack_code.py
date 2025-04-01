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

""" until here just visualization of CellChat outout to understand the data"""



#%%






# %%
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load your data (assuming already loaded as prime_kg_df)
df_subset = prime_kg_casted_df.head(50)  # First 50 rows

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

""" Load the mice knowledge graph """

# Load the saved mouse knowledge graph
mouse_kg_df = pd.read_csv('./data/mouse_translated_prime_kg.csv')

# %%
mouse_kg_df
# %%

""" reminder: i now can use the highly up or down regulated "interesting"
interaction to look relevant stuff up from the new mice kg"""


# --- Highly up/down regulated interactions ---
highly_upregulated_df
#%%

highly_downregulated_df


# %%

# For each ligand-receptor pair in highly_upregulated_df and highly_downregulated_df, you want to: Check if there's a direct or indirect connection in mouse_kg_df.
# Optionally, retrieve a path (e.g., ligand → intermediary → receptor).
# Optionally, annotate with pathway or disease info using public resources.


def find_direct_interactions(df, kg):
    return df[df.apply(lambda row: ((kg['gene1_x'] == row['ligand']) & (kg['gene2_y'] == row['receptor'])).any() |
                                   ((kg['gene1_x'] == row['receptor']) & (kg['gene2_y'] == row['ligand'])).any(), axis=1)]
    
direct_hits_up = find_direct_interactions(highly_upregulated_df, mouse_kg_df)
direct_hits_down = find_direct_interactions(highly_downregulated_df, mouse_kg_df)

# %%
direct_hits_up
# %%
direct_hits_down
# %%
import networkx as nx

# Build the graph
G = nx.from_pandas_edgelist(mouse_kg_df, source='gene1_x', target='gene2_y', edge_attr='relation')

def find_paths(df, G, max_len=2):
    paths = []
    for _, row in df.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']
        if ligand in G and receptor in G:
            try:
                for path in nx.all_simple_paths(G, source=ligand, target=receptor, cutoff=max_len):
                    paths.append({'ligand': ligand, 'receptor': receptor, 'path': path})
            except nx.NetworkXNoPath:
                continue
    return paths

indirect_paths_up = find_paths(highly_upregulated_df, G)
indirect_paths_down = find_paths(highly_downregulated_df, G)

# %%
""" bioservices """
#does not run
from bioservices import KEGG
k = KEGG()

# Get pathways for a gene
def get_kegg_pathways(gene_symbol):
    try:
        res = k.find("genes", gene_symbol)
        if res:
            return res.split("\n")[0]
    except Exception as e:
        return None

# %%
import pandas as pd
import requests

# Loaded data

# Collect unique ligands and receptors
def extract_unique_genes(df):
    return set(df['ligand']).union(set(df['receptor']))

genes_up = extract_unique_genes(highly_upregulated_df)
genes_down = extract_unique_genes(highly_downregulated_df)

# Define enrichment function
def enrichr_query(genes, library='KEGG_2021_Mouse'):
    print(f"Submitting {len(genes)} genes to Enrichr for: {library}")
    ENRICHR_ADD_URL = 'https://maayanlab.cloud/Enrichr/addList'
    ENRICHR_ENRICH_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    
    # Prepare submission
    payload = {
        'list': '\n'.join(genes),
        'description': 'LR pair enrichment'
    }
    response = requests.post(ENRICHR_ADD_URL, files=payload)
    if not response.ok:
        raise Exception('Error submitting gene list to Enrichr')

    user_list_id = response.json()['userListId']

    # Query enrichment results
    params = {'userListId': user_list_id, 'backgroundType': library}
    enrich_response = requests.get(ENRICHR_ENRICH_URL, params=params)
    if not enrich_response.ok:
        raise Exception('Error fetching enrichment results')

    return enrich_response.json()

# Run enrichment for KEGG and Reactome
up_kegg = enrichr_query(genes_up, library='KEGG_2021_Mouse')
up_reactome = enrichr_query(genes_up, library='Reactome_2016')

down_kegg = enrichr_query(genes_down, library='KEGG_2021_Mouse')
down_reactome = enrichr_query(genes_down, library='Reactome_2016')

# Helper to parse results
def display_top_terms(results, top_n=10):
    df = pd.DataFrame(results['KEGG_2021_Mouse'] if 'KEGG_2021_Mouse' in results else list(results.values())[0],
                      columns=['Rank', 'Term', 'P-value', 'Z-score', 'Combined score', 'Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'])
    return df.head(top_n)

# View results
# View results
print("\n🔺 Top KEGG Pathways (UP):")
print(display_top_terms(up_kegg))

print("\n🔻 Top KEGG Pathways (DOWN):")
print(display_top_terms(down_kegg))

print("\n🔺 Top Reactome Pathways (UP):")
print(display_top_terms(up_reactome))

print("\n🔻 Top Reactome Pathways (DOWN):")
print(display_top_terms(down_reactome))

# %%

#🔺 Top KEGG Pathways (Upregulated LR Pairs):


#🔻 Top KEGG Pathways (Downregulated LR Pairs):


#%%

#Collecting ligands and receptors

#Querying Enrichr

#Parsing results

#Visualizing top pathways

""" reminder: i now can use the highly up or down regulated "interesting"
interaction to look relevant stuff up from the new mice kg"""


# --- Highly up/down regulated interactions ---
highly_upregulated_df


highly_downregulated_df



# %%
import matplotlib.pyplot as plt
import pandas as pd

# Sample data from previous enrichment (mocked here for reproducibility)
# You would replace this with the actual data from the enrichment results
top_up_kegg_df = pd.DataFrame({
    'Term': ['Antigen processing', 'Cell adhesion', 'Allograft rejection', 'Graft vs host', 'Type I diabetes',
             'Autoimmune thyroid', 'Viral myocarditis', 'T-cell leukemia', 'Epstein-Barr', 'Phagosome'],
    'Combined score': [6641.97, 2043.48, 4318.77, 4318.77, 2213.47, 1777.17, 1542.55, 1181.00, 965.46, 843.00]
})

top_down_kegg_df = pd.DataFrame({
    'Term': ['ECM-receptor interaction', 'Viral interaction', 'Chemokine signaling', 'Cytokine-cytokine',
             'Toxoplasmosis', 'Hematopoietic lineage', 'Focal adhesion', 'Rheumatoid arthritis',
             'Antigen processing', 'Cell adhesion'],
    'Combined score': [5346.10, 2849.52, 1046.25, 745.43, 888.33, 745.83, 373.18, 524.96, 508.33, 310.85]
})

# Determine the maximum combined score for consistent Y-axis
max_score = max(top_up_kegg_df['Combined score'].max(), top_down_kegg_df['Combined score'].max())

# Plot Upregulated Pathways
plt.figure(figsize=(10, 6))
df_sorted_up = top_up_kegg_df.sort_values('Combined score', ascending=True)
plt.barh(df_sorted_up['Term'], df_sorted_up['Combined score'], edgecolor='black', color='red')
plt.xlabel("Combined score")
plt.title("🔺 Top KEGG Pathways (Upregulated LR Pairs)")
plt.xlim([0, max_score])
plt.tight_layout()
plt.show()

# Plot Downregulated Pathways
plt.figure(figsize=(10, 6))
df_sorted_down = top_down_kegg_df.sort_values('Combined score', ascending=True)
plt.barh(df_sorted_down['Term'], df_sorted_down['Combined score'], edgecolor='black', color='skyblue')
plt.xlabel("Combined score")
plt.title("🔻 Top KEGG Pathways (Downregulated LR Pairs)")
plt.xlim([0, max_score])
plt.tight_layout()
plt.show()

# %%
