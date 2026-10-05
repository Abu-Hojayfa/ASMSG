import os
import json
import numpy as np
import pandas as pd

def estimate_pmid_year(pmid_str):
    """
    Estimates publication year from PubMed ID using empirical PMID-to-year mapping thresholds.
    """
    try:
        pmid = int(str(pmid_str).strip())
    except (ValueError, TypeError):
        return 2015 # default median fallback
        
    if pmid < 10000000:
        return 1995
    elif pmid < 15000000:
        return 2004
    elif pmid < 20000000:
        return 2009
    elif pmid < 25000000:
        return 2014
    elif pmid < 30000000:
        return 2018
    elif pmid < 33000000:
        return 2020
    elif pmid < 36000000:
        return 2022
    else:
        return 2024

def build_temporal_split(cutoff_year=2020):
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    edges_csv = os.path.join(base_dir, 'datasets', 'asmsg_clean', 'asmsg_edges.csv')
    out_dir = os.path.join(base_dir, 'datasets', 'asmsg_clean')
    
    if not os.path.exists(edges_csv):
        raise FileNotFoundError(f"Edges table not found at {edges_csv}")
        
    print(f"\n--- Generating Temporal Validation Split (Cutoff Year T = {cutoff_year}) ---")
    df = pd.read_csv(edges_csv)
    
    train_indices = []
    test_indices = []
    
    for idx, row in df.iterrows():
        pmid_list = str(row['PMID_List']).split(';')
        years = [estimate_pmid_year(p) for p in pmid_list if p.isdigit()]
        
        if not years:
            min_year = 2015
        else:
            min_year = min(years)
            
        if min_year < cutoff_year:
            train_indices.append(idx)
        else:
            test_indices.append(idx)
            
    print(f"Total Edges: {len(df)}")
    print(f" - Historical Training Edges (< {cutoff_year}): {len(train_indices)} ({len(train_indices)/len(df)*100:.1f}%)")
    print(f" - Prospective Test Edges (>= {cutoff_year}): {len(test_indices)} ({len(test_indices)/len(df)*100:.1f}%)")
    
    temporal_data = {
        'cutoff_year': cutoff_year,
        'train_edge_count': len(train_indices),
        'test_edge_count': len(test_indices),
        'train_edge_indices': train_indices,
        'test_edge_indices': test_indices
    }
    
    out_file = os.path.join(out_dir, f'asmsg_temporal_split_cutoff_{cutoff_year}.json')
    with open(out_file, 'w', encoding='utf-8') as f:
        json.dump(temporal_data, f, indent=2)
        
    print(f"Saved temporal split to {out_file}")
    return temporal_data

if __name__ == "__main__":
    build_temporal_split(cutoff_year=2020)
