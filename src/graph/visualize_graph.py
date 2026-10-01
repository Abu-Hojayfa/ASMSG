import torch
import networkx as nx
import matplotlib.pyplot as plt
import os

print("Loading Graph...")
graph = torch.load('d:/fydp/features/asmsg_hetero_graph.pt', weights_only=False)

G = nx.Graph()

# We will sample a very small subset so the visualization is readable
# Let's take the first 50 RNA-Disease edges and 50 RNA-Drug edges
rna_disease_edges = graph['rna', 'associates_with', 'disease'].edge_index
rna_drug_edges = graph['rna', 'targets', 'drug'].edge_index

num_samples = 40

print("Building NetworkX Graph...")
# Add RNA-Disease edges
for i in range(num_samples):
    r_id = f"RNA_{rna_disease_edges[0, i].item()}"
    d_id = f"Disease_{rna_disease_edges[1, i].item()}"
    G.add_node(r_id, type='RNA')
    G.add_node(d_id, type='Disease')
    G.add_edge(r_id, d_id, type='associates_with')

# Add RNA-Drug edges
for i in range(num_samples):
    r_id = f"RNA_{rna_drug_edges[0, i].item()}"
    dr_id = f"Drug_{rna_drug_edges[1, i].item()}"
    G.add_node(r_id, type='RNA')
    G.add_node(dr_id, type='Drug')
    G.add_edge(r_id, dr_id, type='targets')

# Color map
color_map = []
for node, data in G.nodes(data=True):
    if data['type'] == 'RNA':
        color_map.append('#ff9999') # Red
    elif data['type'] == 'Disease':
        color_map.append('#99ccff') # Blue
    else:
        color_map.append('#99ff99') # Green

print("Plotting graph...")
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.5, iterations=50)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=500, edgecolors='black')

# Draw edges
nx.draw_networkx_edges(G, pos, alpha=0.5, width=2.0)

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=8, font_family='sans-serif')

plt.title("ASMSG Heterogeneous Graph Sub-sample\n(Red=RNA, Blue=Disease, Green=Drug)", fontsize=16)
plt.axis('off')

# Save to the specific artifact directory so the user can see it in chat!
out_path = r"C:\Users\abuho\.gemini\antigravity-ide\brain\d9353bae-82ea-42de-9f66-d27771628a01\graph_visualization.png"
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"Saved visualization to {out_path}")
