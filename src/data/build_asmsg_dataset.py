import os
import json
import numpy as np
import pandas as pd

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
    # REVISED DISEASE ONTOLOGY & ENTITY RESOLUTION ENGINE
    # ----------------─────────────────────────────────────────────────────
    print("\nBuilding dataset-wide Disease Ontology inheritance maps...")
    
    # Build dataset-wide cross-row mapping tables
    name_to_doid = {}
    name_to_mesh = {}
    doid_to_canonical_name = {}
    mesh_to_canonical_name = {}
    
    for _, row in df.dropna(subset=['DO ID']).iterrows():
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        doid = str(row['DO ID']).strip()
        if doid and doid.lower() not in ['nan', 'none']:
            norm_doid = doid if doid.startswith('DOID:') else f"DOID:{doid}"
            name_to_doid[dname_lower] = norm_doid
            if norm_doid not in doid_to_canonical_name:
                doid_to_canonical_name[norm_doid] = dname

    for _, row in df.dropna(subset=['MeSH ID']).iterrows():
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        mesh = str(row['MeSH ID']).strip()
        if mesh and mesh.lower() not in ['nan', 'none']:
            norm_mesh = mesh if mesh.startswith('MESH:') else f"MESH:{mesh}"
            name_to_mesh[dname_lower] = norm_mesh
            if norm_mesh not in mesh_to_canonical_name:
                mesh_to_canonical_name[norm_mesh] = dname

    print(f" - Mapped {len(name_to_doid)} disease names to official DO IDs.")
    print(f" - Mapped {len(name_to_mesh)} disease names to official MeSH IDs.")
    
    def resolve_disease(row):
        dname = str(row['Disease Name']).strip()
        dname_lower = dname.lower()
        
        # 1. Primary: DO ID (explicit or dataset-inherited)
        doid = str(row['DO ID']).strip() if pd.notnull(row['DO ID']) else ""
        if not doid or doid.lower() in ['nan', 'none']:
            doid = name_to_doid.get(dname_lower, "")
        else:
            doid = doid if doid.startswith('DOID:') else f"DOID:{doid}"
            
        if doid and doid.lower() not in ['nan', 'none']:
            canonical_name = doid_to_canonical_name.get(doid, dname)
            return doid, canonical_name, doid.replace('DOID:', ''), name_to_mesh.get(dname_lower, "")
            
        # 2. Secondary: MeSH ID (explicit or dataset-inherited)
        mesh = str(row['MeSH ID']).strip() if pd.notnull(row['MeSH ID']) else ""
        if not mesh or mesh.lower() in ['nan', 'none']:
            mesh = name_to_mesh.get(dname_lower, "")
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
    # OUTPUT MASTER CSV TABLES
    # ----------------─────────────────────────────────────────────────────
    
    # 1. Build Nodes CSV (ncRNA Master Table)
    print("\nBuilding asmsg_nodes.csv...")
    node_data = df_clean.drop_duplicates(subset=['RNA Symbol'])[['RNA Symbol', 'RNA Type']].copy()
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
    
    # 3. Build Edges CSV with PMID Evidence Weighting
    print("Building asmsg_edges.csv with continuous PMID evidence weighting...")
    
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
    
    # Compute continuous log-scaled edge weights: w_ij = log(1 + E_ij) / max(log(1 + E_ij))
    max_log_ev = np.log1p(edge_df['Evidence_Count'].max())
    edge_df['Edge_Weight'] = (np.log1p(edge_df['Evidence_Count']) / max_log_ev).round(4)
    
    edge_file = os.path.join(out_dir, 'asmsg_edges.csv')
    edge_df.to_csv(edge_file, index=False)
    
    print("\n" + "="*50)
    print("ASMSG REVISED DATASET CONSTRUCTION COMPLETE!")
    print(f"Unique ncRNA Nodes: {len(node_data)}")
    print(f"Unique Standardized Disease Nodes: {len(disease_group)}")
    print(f" - DO ID Nodes: {len(disease_group[disease_group['Disease_ID'].str.startswith('DOID')])}")
    print(f" - MeSH Nodes: {len(disease_group[disease_group['Disease_ID'].str.startswith('MESH')])}")
    print(f" - Isolated Name Nodes: {len(disease_group[disease_group['Disease_ID'].str.startswith('NAME')])}")
    print(f"Unique ncRNA-Disease Edges: {len(edge_df)}")
    print(f"Max Evidence Count for single edge: {edge_df['Evidence_Count'].max()}")
    print(f"Saved to: {out_dir}")
    print("="*50)

if __name__ == "__main__":
    build_clean_dataset()
