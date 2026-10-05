import os
import re
import json
import itertools
import numpy as np
import pandas as pd

def compute_head_tail_clamping(seq, max_len=1022):
    """
    Clamps long RNA sequences while preserving both 5'-end seed region and 3'-end regulatory domain.
    If len(seq) > 1022, takes first 511 nt + last 511 nt.
    """
    seq_str = str(seq).upper().replace('T', 'U')
    if len(seq_str) <= max_len:
        return seq_str
    half = max_len // 2
    return seq_str[:half] + seq_str[-half:]

def compute_3mer_similarity_edges(node_data, threshold=0.85, top_k=5):
    """
    Computes sequence similarity edges between ncRNAs based on 3-mer frequency cosine similarity
    to connect degree-1 leaf nodes and enrich graph topology.
    """
    print("\nComputing 3-mer sequence similarity edges among ncRNAs...")
    nucleotides = ['A', 'C', 'G', 'U']
    kmers = [''.join(p) for p in itertools.product(nucleotides, repeat=3)]
    kmer_to_idx = {k: i for i, k in enumerate(kmers)}
    
    sequences = node_data['Sequence'].tolist()
    symbols = node_data['RNA Symbol'].tolist()
    matrix = np.zeros((len(sequences), len(kmers)), dtype=np.float32)
    
    for row_idx, seq in enumerate(sequences):
        seq_clean = seq.upper().replace('T', 'U')
        seq_len = len(seq_clean)
        if seq_len < 3:
            continue
        counts = {}
        for i in range(seq_len - 2):
            kmer = seq_clean[i : i + 3]
            if kmer in kmer_to_idx:
                counts[kmer] = counts.get(kmer, 0) + 1
        total = sum(counts.values())
        if total > 0:
            for kmer, cnt in counts.items():
                matrix[row_idx, kmer_to_idx[kmer]] = cnt / total
                
    # Normalize rows for cosine similarity
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    norm_matrix = matrix / norms
    
    sim_edges = []
    chunk_size = 2000
    for i in range(0, len(symbols), chunk_size):
        chunk = norm_matrix[i : i + chunk_size]
        sim_chunk = np.dot(chunk, norm_matrix.T)
        
        for r_local in range(len(chunk)):
            r_global = i + r_local
            sim_scores = sim_chunk[r_local]
            sim_scores[r_global] = 0.0 # Zero self-similarity
            
            top_indices = np.argsort(sim_scores)[-top_k:]
            for target_idx in top_indices:
                score = sim_scores[target_idx]
                if score >= threshold:
                    sim_edges.append({
                        'RNA_Symbol_1': symbols[r_global],
                        'RNA_Symbol_2': symbols[target_idx],
                        'Similarity_Score': round(float(score), 4)
                    })
                    
    sim_df = pd.DataFrame(sim_edges).drop_duplicates()
    print(f" - Generated {len(sim_df)} ncRNA-ncRNA sequence similarity edges (threshold >= {threshold}).")
    return sim_df

def build_clean_dataset():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    excel_path = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'alldata.xlsx')
    fasta_path = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'matched_sequences.fa')
    out_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean')
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
            
    print("Loading matched sequences...")
    matched_seqs = {}
    current_name = None
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_name = line[1:]
            elif current_name:
                matched_seqs[current_name] = line
                current_name = None
                
    valid_names = set(matched_seqs.keys())
    print(f"Total verified RNA nodes with sequences: {len(valid_names)}")
    
    print("Loading RNADisease database (alldata.xlsx)...")
    df = pd.read_excel(excel_path)
    
    # Filter for Human and ncRNA
    df = df[df['specise'].astype(str).str.contains("Homo sapiens", na=False, case=False)]
    df = df[~df['RNA Type'].astype(str).str.contains("mRNA", na=False, case=False)]
    
    print(f"Total human ncRNA rows: {len(df)}")
    
    # Keep ONLY rows where the RNA Symbol is in our valid, matched list
    df_clean = df[df['RNA Symbol'].isin(valid_names)].copy()
    print(f"Human ncRNA rows after sequence verification: {len(df_clean)}")
    
    # ----------------─────────────────────────────────────────────────────
    # STRICT ONTOLOGY DISEASE RESOLUTION (DO ID & MeSH ONLY — OPTION B)
    # ----------------─────────────────────────────────────────────────────
    print("\nBuilding dataset-wide Disease Ontology & Synonym inheritance maps...")
    
    name_to_doid = {}
    name_to_mesh = {}
    doid_to_canonical_name = {}
    mesh_to_canonical_name = {}
    
    def clean_name(s):
        s_clean = str(s).strip().lower()
        s_clean = re.sub(r'\s+', ' ', s_clean)
        s_clean = re.sub(r' (disease|syndrome|carcinoma|cancer|neoplasm)$', '', s_clean)
        return s_clean
        
    for _, row in df.dropna(subset=['DO ID']).iterrows():
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        dname_c = clean_name(dname)
        doid = str(row['DO ID']).strip()
        if doid and doid.lower() not in ['nan', 'none']:
            norm_doid = doid if doid.startswith('DOID:') else f"DOID:{doid}"
            name_to_doid[dname_lower] = norm_doid
            name_to_doid[dname_c] = norm_doid
            if norm_doid not in doid_to_canonical_name:
                doid_to_canonical_name[norm_doid] = dname

    for _, row in df.dropna(subset=['MeSH ID']).iterrows():
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        dname_c = clean_name(dname)
        mesh = str(row['MeSH ID']).strip()
        if mesh and mesh.lower() not in ['nan', 'none']:
            norm_mesh = mesh if mesh.startswith('MESH:') else f"MESH:{mesh}"
            name_to_mesh[dname_lower] = norm_mesh
            name_to_mesh[dname_c] = norm_mesh
            if norm_mesh not in mesh_to_canonical_name:
                mesh_to_canonical_name[norm_mesh] = dname

    def resolve_disease_strict(row):
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        dname_c = clean_name(dname)
        
        # 1. Primary: DO ID (explicit or dataset-inherited)
        doid = str(row['DO ID']).strip() if pd.notnull(row['DO ID']) else ""
        if not doid or doid.lower() in ['nan', 'none']:
            doid = name_to_doid.get(dname_lower, name_to_doid.get(dname_c, ""))
        else:
            doid = doid if doid.startswith('DOID:') else f"DOID:{doid}"
            
        if doid and doid.lower() not in ['nan', 'none']:
            canonical_name = doid_to_canonical_name.get(doid, dname)
            return doid, canonical_name, doid.replace('DOID:', ''), name_to_mesh.get(dname_lower, "")
            
        # 2. Secondary: MeSH ID (explicit or dataset-inherited)
        mesh = str(row['MeSH ID']).strip() if pd.notnull(row['MeSH ID']) else ""
        if not mesh or mesh.lower() in ['nan', 'none']:
            mesh = name_to_mesh.get(dname_lower, name_to_mesh.get(dname_c, ""))
        else:
            mesh = mesh if mesh.startswith('MESH:') else f"MESH:{mesh}"
            
        if mesh and mesh.lower() not in ['nan', 'none']:
            canonical_name = mesh_to_canonical_name.get(mesh, dname)
            return mesh, canonical_name, "", mesh.replace('MESH:', '')
            
        # Filter out unmapped strings strictly for Option B
        return None, None, None, None

    resolved_records = df_clean.apply(resolve_disease_strict, axis=1)
    df_clean['Disease_ID'] = [r[0] for r in resolved_records]
    df_clean['Canonical_Disease_Name'] = [r[1] for r in resolved_records]
    df_clean['Resolved_DO_ID'] = [r[2] for r in resolved_records]
    df_clean['Resolved_MeSH_ID'] = [r[3] for r in resolved_records]
    
    # Strictly remove unmapped string rows
    df_clean = df_clean.dropna(subset=['Disease_ID']).copy()
    print(f"Association rows after 100% strict ontology disease resolution: {len(df_clean)}")
    
    # --------------------------------─────────────────────────────────────
    # BIOTYPE STANDARDIZATION
    # ----------------─────────────────────────────────────────────────────
    def standardize_biotype(btype):
        b = str(btype).strip().lower()
        if 'mirna' in b:
            return 'miRNA'
        if 'lncrna' in b:
            return 'lncRNA'
        if 'circrna' in b:
            return 'circRNA'
        if 'pirna' in b:
            return 'piRNA'
        if 'snorna' in b:
            return 'snoRNA'
        return 'other_ncRNA'
        
    df_clean['Standardized_RNA_Type'] = df_clean['RNA Type'].apply(standardize_biotype)
    
    # ----------------─────────────────────────────────────────────────────
    # OUTPUT SINGLE PERFECTED MASTER CSV TABLES (asmsg_clean/)
    # ----------------─────────────────────────────────────────────────────
    
    # 1. Build Nodes CSV (ncRNA Master Table)
    print("\nBuilding asmsg_nodes.csv with Head-Tail Dual Window Clamping...")
    node_data = df_clean.drop_duplicates(subset=['RNA Symbol'])[['RNA Symbol', 'Standardized_RNA_Type']].copy()
    node_data.rename(columns={'Standardized_RNA_Type': 'RNA Type'}, inplace=True)
    
    # Apply Head-Tail dual-window sequence clamping
    raw_seqs = [matched_seqs[sym] for sym in node_data['RNA Symbol']]
    clamped_seqs = [compute_head_tail_clamping(seq, max_len=1022) for seq in raw_seqs]
    
    node_data['Sequence'] = clamped_seqs
    node_data['Raw_Seq_Length'] = [len(s) for s in raw_seqs]
    node_data['Clamped_Seq_Length'] = [len(s) for s in clamped_seqs]
    node_data['Is_Clamped'] = node_data['Raw_Seq_Length'] > 1022
    
    node_file = os.path.join(out_dir, 'asmsg_nodes.csv')
    node_data.to_csv(node_file, index=False)
    
    # 2. Build Diseases CSV (Disease Master Table — 100% DOID/MeSH Curated)
    print("Building asmsg_diseases.csv (100% DOID/MeSH Curated)...")
    disease_group = df_clean.groupby('Disease_ID').agg(
        Disease_Name=('Canonical_Disease_Name', 'first'),
        DO_ID=('Resolved_DO_ID', 'first'),
        MeSH_ID=('Resolved_MeSH_ID', 'first'),
        Total_Edges=('RNA Symbol', 'count')
    ).reset_index()
    disease_file = os.path.join(out_dir, 'asmsg_diseases.csv')
    disease_group.to_csv(disease_file, index=False)
    
    # 3. Build Edges CSV with Dual Continuous Evidence Weighting
    print("Building asmsg_edges.csv with composite dual weighting...")
    
    def process_pmids(series):
        pmids = set()
        for val in series.dropna():
            s_val = str(val).strip()
            if s_val and s_val != 'nan':
                for p in s_val.replace(';', ',').split(','):
                    p_clean = p.strip().split('.')[0]
                    if p_clean.isdigit():
                        pmids.add(p_clean)
        return list(pmids)

    def process_scores(series):
        valid_scores = []
        for val in series.dropna():
            try:
                s_float = float(val)
                if not np.isnan(s_float):
                    valid_scores.append(s_float)
            except (ValueError, TypeError):
                pass
        return round(float(max(valid_scores)), 4) if valid_scores else 0.5000

    edge_summary = []
    grouped_edges = df_clean.groupby(['RNA Symbol', 'Disease_ID'])
    
    for (rna_sym, dis_id), group in grouped_edges:
        pmids = process_pmids(group['PMID'])
        ev_count = max(len(pmids), len(group))
        pmid_str = ";".join(pmids) if pmids else "UNSPECIFIED"
        conf_score = process_scores(group['score'])
        
        edge_summary.append({
            'RNA Symbol': rna_sym,
            'Disease_ID': dis_id,
            'Evidence_Count': ev_count,
            'PMID_List': pmid_str,
            'Confidence_Score': conf_score
        })
        
    edge_df = pd.DataFrame(edge_summary)
    
    # Dual continuous weight formula: w_ij = 0.5 * w_PMID + 0.5 * Confidence_Score
    max_log_ev = np.log1p(edge_df['Evidence_Count'].max())
    pmid_weight = np.log1p(edge_df['Evidence_Count']) / max_log_ev
    edge_df['Edge_Weight'] = (0.5 * pmid_weight + 0.5 * edge_df['Confidence_Score']).round(4)
    
    edge_file = os.path.join(out_dir, 'asmsg_edges.csv')
    edge_df.to_csv(edge_file, index=False)
    
    # 4. Compute ncRNA Topological Sequence Similarity Edges
    sim_edges_df = compute_3mer_similarity_edges(node_data, threshold=0.85, top_k=5)
    sim_edges_file = os.path.join(out_dir, 'ncrna_sequence_similarity_edges.csv')
    sim_edges_df.to_csv(sim_edges_file, index=False)
    
    print("\n" + "="*50)
    print("ASMSG SINGLE MASTER DATASET (asmsg_clean/) OPTION B COMPLETE!")
    print(f"ncRNA Nodes: {len(node_data)} (Clamped at 1022 nt via Head-Tail Dual Window: {node_data['Is_Clamped'].sum()})")
    print(f"Standardized Disease Nodes: {len(disease_group)} (DOID: {len(disease_group[disease_group['Disease_ID'].str.startswith('DOID:')])}, MeSH: {len(disease_group[disease_group['Disease_ID'].str.startswith('MESH:')])}, Raw Strings: 0)")
    print(f"ncRNA-Disease Associations: {len(edge_df)}")
    print(f"ncRNA-ncRNA Sequence Similarity Edges: {len(sim_edges_df)}")
    print(f"Edge Weight Mean: {edge_df['Edge_Weight'].mean():.4f}, Std: {edge_df['Edge_Weight'].std():.4f}, Min: {edge_df['Edge_Weight'].min():.4f}, Max: {edge_df['Edge_Weight'].max():.4f}")
    print("="*50)

if __name__ == "__main__":
    build_clean_dataset()
