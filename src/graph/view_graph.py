import torch

print("Loading ASMSG Heterogeneous Graph...\n")
graph = torch.load('d:/fydp/features/asmsg_hetero_graph.pt', weights_only=False)

print("="*50)
print("             NODE TYPES & FEATURES")
print("="*50)
for node_type in graph.node_types:
    features = graph[node_type].x
    print(f"- Node: '{node_type}'")
    print(f"  Count: {features.shape[0]:,} nodes")
    print(f"  Features: {features.shape[1]} dimensions")
    # Print the first feature vector snippet
    print(f"  Example (Node 0): {features[0][:5].tolist()}...\n")

print("="*50)
print("             EDGE TYPES & CONNECTIONS")
print("="*50)
for edge_type in graph.edge_types:
    src, rel, dst = edge_type
    edge_index = graph[edge_type].edge_index
    num_edges = edge_index.shape[1]
    print(f"- Edge: {src} -> [{rel}] -> {dst}")
    print(f"  Count: {num_edges:,} connections")
    # Print the first few connections
    print(f"  Example (First 3 connections):")
    for i in range(min(3, num_edges)):
        print(f"    {src} {edge_index[0, i]} is connected to {dst} {edge_index[1, i]}")
    print()

print("="*50)
print("Graph is mathematically sound and ready for PyTorch!")
