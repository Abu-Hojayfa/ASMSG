import os
import re
import json
import numpy as np
import pandas as pd

def build_clean_dataset():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    excel_path = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'alldata.xlsx')
    fasta_path = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'matched_sequences.fa')
    out_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean')
    core_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean_core')
    
    for d in [out_dir, core_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
            
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
    # ADVANCED DISEASE ONTOLOGY & SYNONYM RESOLUTION ENGINE
    # ----------------─────────────────────────────────────────────────────
    print("\nBuilding dataset-wide Disease Ontology & Synonym inheritance maps...")
    
    name_to_doid = {}
    name_to_mesh = {}
    doid_to_canonical_name = {}
    mesh_to_canonical_name = {}
    
    def clean_name(s):
        s_clean = str(s).strip().lower()
        s_clean = re.sub(r'\s+', ' ', s_clean)
        # Remove common trailing noise words for matching
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

    print(f" - Mapped {len(name_to_doid)} disease names/synonyms to official DO IDs.")
    print(f" - Mapped {len(name_to_mesh)} disease names/synonyms to official MeSH IDs.")
    
    def resolve_disease(row):
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
            
        # 3. Tertiary: Normalized Phenotype String
        return f"NAME:{dname_lower}", dname, "", ""

    resolved_records = df_clean.apply(resolve_disease, axis=1)
    df_clean['Disease_ID'] = [r[0] for r in resolved_records]
    df_clean['Canonical_Disease_Name'] = [r[1] for r in resolved_records]
    df_clean['Resolved_DO_ID'] = [r[2] for r in resolved_records]
    df_clean['Resolved_MeSH_ID'] = [r[3] for r in resolved_records]
    
    # ----------------─────────────────────────────────────────────────────
    # BIOTYPE STANDARDIZATION
    # ----------------─────────────────────────────────────────────────────
    def standardize_biotype(btype):
        b = str(btype).strip()
        b_lower = b.lower()
        if 'mirna' in b_lower:
            return 'miRNA'
        if 'lncrna' in b_lower:
            return 'lncRNA'
        if 'circrna' in b_lower:
            return 'circRNA'
        if 'pirna' in b_lower:
            return 'piRNA'
        if 'snorna' in b_lower:
            return 'snoRNA'
        return 'other_ncRNA'
        
    df_clean['Standardized_RNA_Type'] = df_clean['RNA Type'].apply(standardize_biotype)
    
    # ----------------─────────────────────────────────────────────────────
    # OUTPUT MASTER CSV TABLES (asmsg_clean/)
    # ----------------─────────────────────────────────────────────────────
    
    # 1. Build Nodes CSV (ncRNA Master Table)
    print("\nBuilding asmsg_nodes.csv...")
    node_data = df_clean.drop_duplicates(subset=['RNA Symbol'])[['RNA Symbol', 'Standardized_RNA_Type']].copy()
    node_data.rename(columns={'Standardized_RNA_Type': 'RNA Type'}, inplace=True)
    node_data['Sequence'] = node_data['RNA Symbol'].map(matched_seqs)
    node_data['Seq_Length'] = node_data['Sequence'].str.len()
    node_file = os.path.join(out_dir, 'asmsg_nodes.csv')
    node_data.to_csv(node_file, index=False)
    
    # 2. Build Diseases CSV (Disease Master Table)
    print("Building asmsg_diseases.csv...")
    disease_group = df_clean.groupby('Disease_ID').agg(
        Disease_Name=('Canonical_Disease_Name', 'first'),
        DO_ID=('Resolved_DO_ID', 'first'),
        MeSH_ID=('Resolved_MeSH_ID', 'first'),
        Total_Edges=('RNA Symbol', 'count')
    ).reset_index()
    disease_file = os.path.join(out_dir, 'asmsg_diseases.csv')
    disease_group.to_csv(disease_file, index=False)
    
    # 3. Build Edges CSV with Dual Continuous Evidence Weighting
    print("Building asmsg_edges.csv with dual PMID & confidence score weighting...")
    
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
    
    # ----------------─────────────────────────────────────────────────────
    # BUILD HIGH-CONFIDENCE CORE GRAPH VARIANT (asmsg_clean_core/)
    # ----------------─────────────────────────────────────────────────────
    print("\nBuilding high-confidence core dataset variant (asmsg_clean_core/)...")
    
    # Core filter: Exclude unmapped NAME: disease nodes
    core_disease_ids = set(disease_group[~disease_group['Disease_ID'].str.startswith('NAME:')]['Disease_ID'])
    core_edges_df = edge_df[edge_df['Disease_ID'].isin(core_disease_ids)].copy()
    
    core_rna_symbols = set(core_edges_df['RNA Symbol'].unique())
    core_nodes_df = node_data[node_data['RNA Symbol'].isin(core_rna_symbols)].copy()
    core_diseases_df = disease_group[disease_group['Disease_ID'].isin(core_disease_ids)].copy()
    
    core_nodes_df.to_csv(os.path.join(core_dir, 'asmsg_nodes.csv'), index=False)
    core_diseases_df.to_csv(os.path.join(core_dir, 'asmsg_diseases.csv'), index=False)
    core_edges_df.to_csv(os.path.join(core_dir, 'asmsg_edges.csv'), index=False)
    
    print("\n" + "="*50)
    print("ASMSG REVISED DATASET RE-CONSTRUCTION COMPLETE!")
    print(f"FULL GRAPH (asmsg_clean/):")
    print(f" - ncRNA Nodes: {len(node_data)}")
    print(f" - Disease Nodes: {len(disease_group)} (DOID: {len(disease_group[disease_group['Disease_ID'].str.startswith('DOID:')])}, MeSH: {len(disease_group[disease_group['Disease_ID'].str.startswith('MESH:')])}, NAME: {len(disease_group[disease_group['Disease_ID'].str.startswith('NAME:')])})")
    print(f" - Edges: {len(edge_df)}")
    print(f" - Edge Weights Mean: {edge_df['Edge_Weight'].mean():.4f}, Std: {edge_df['Edge_Weight'].std():.4f}, Min: {edge_df['Edge_Weight'].min():.4f}, Max: {edge_df['Edge_Weight'].max():.4f}")
    
    print(f"\nCORE GRAPH (asmsg_clean_core/):")
    print(f" - ncRNA Nodes: {len(core_nodes_df)}")
    print(f" - Standardized Disease Nodes: {len(core_diseases_df)}")
    print(f" - High-Confidence Edges: {len(core_edges_df)}")
    print("="*50)

if __name__ == "__main__":
    build_clean_dataset()
