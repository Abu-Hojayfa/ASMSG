import os
import re
import json
import gzip
import pandas as pd
from collections import defaultdict

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

def parse_fasta_file(file_path, db_name):
    """
    Parses a FASTA or FASTA.GZ file into a dictionary of {key: sequence}.
    """
    seq_dict = {}
    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} not found. Skipping.")
        return seq_dict
        
    print(f"Parsing FASTA database [{db_name}] from {os.path.basename(file_path)}...")
    
    is_gz = file_path.endswith('.gz')
    open_fn = lambda p: gzip.open(p, 'rt', encoding='utf-8', errors='ignore') if is_gz else open(p, 'r', encoding='utf-8', errors='ignore')
    
    with open_fn(file_path) as f:
        current_keys = []
        current_seq_lines = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_keys and current_seq_lines:
                    full_seq = "".join(current_seq_lines).upper().replace('T', 'U')
                    # ENFORCE MINIMUM SEQUENCE LENGTH GUARD (>= 15 nt) to filter single/triple nt garbage
                    if len(full_seq) >= 15:
                        clamped_seq = compute_head_tail_clamping(full_seq, max_len=1022)
                        for k in current_keys:
                            if k and k not in seq_dict:
                                seq_dict[k] = clamped_seq
                                
                current_keys = []
                current_seq_lines = []
                header = line[1:].strip()
                
                # Extract lookup keys from header based on database format
                if db_name == 'mirbase_mature':
                    # Header format: >hsa-let-7a-5p MIMAT0000062 ...
                    parts = header.split()
                    current_keys.append(parts[0].lower())
                    if len(parts) > 1:
                        current_keys.append(parts[1].lower())
                elif db_name == 'mirbase_hairpin':
                    # Header format: >hsa-let-7a-1 MI0000060 ...
                    # Prefix with 'hairpin:' to NEVER overwrite mature miRNAs!
                    parts = header.split()
                    current_keys.append(f"hairpin:{parts[0].lower()}")
                    if len(parts) > 1:
                        current_keys.append(f"hairpin:{parts[1].lower()}")
                elif db_name in ['gencode', 'ensembl']:
                    parts = header.split('|')
                    for p in parts:
                        p_clean = p.strip()
                        if p_clean:
                            current_keys.append(p_clean.lower())
                            # Strip version numbers e.g. ENSG00000228630.1 -> ENSG00000228630
                            if '.' in p_clean:
                                current_keys.append(p_clean.split('.')[0].lower())
                else:
                    parts = header.split()
                    current_keys.append(parts[0].lower())
            else:
                current_seq_lines.append(line)
                
        if current_keys and current_seq_lines:
            full_seq = "".join(current_seq_lines).upper().replace('T', 'U')
            if len(full_seq) >= 15:
                clamped_seq = compute_head_tail_clamping(full_seq, max_len=1022)
                for k in current_keys:
                    if k and k not in seq_dict:
                        seq_dict[k] = clamped_seq
                        
    print(f" - Loaded {len(seq_dict)} sequence keys from [{db_name}].")
    return seq_dict

def match_sequences():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    excel_path = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'alldata.xlsx')
    out_fasta = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'matched_sequences.fa')
    out_json = os.path.join(base_dir, 'datasets', 'rnadisease_v4', 'unique_ncrnas.json')
    
    print("Loading RNADisease v4.0 target ncRNA names...")
    df = pd.read_excel(excel_path)
    df = df[df['specise'].astype(str).str.contains("Homo sapiens", na=False, case=False)]
    df = df[~df['RNA Type'].astype(str).str.contains("mRNA", na=False, case=False)]
    
    target_names = sorted(list(set(df['RNA Symbol'].dropna().astype(str).str.strip())))
    print(f"Total target human ncRNA symbols: {len(target_names)}")
    
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(target_names, f, indent=2)
        
    print("\nParsing reference sequence databases with priority hierarchy...")
    
    # 1. Parse mature miRNAs FIRST (Highest priority for miRNAs)
    mature_path = os.path.join(base_dir, 'datasets', 'mirbase_sequences', 'mature.fa')
    mature_seqs = parse_fasta_file(mature_path, 'mirbase_mature')
    
    # 2. Parse GENCODE lncRNAs
    gencode_path = os.path.join(base_dir, 'datasets', 'ncrna_sequences', 'gencode_lncrna.fa.gz')
    gencode_seqs = parse_fasta_file(gencode_path, 'gencode')
    
    # 3. Parse Ensembl ncRNAs
    ensembl_path = os.path.join(base_dir, 'datasets', 'ncrna_sequences', 'ensembl_ncrna.fa.gz')
    ensembl_seqs = parse_fasta_file(ensembl_path, 'ensembl')
    
    # 4. Parse hairpin precursor miRNAs with explicit 'hairpin:' namespace guard
    hairpin_path_gz = os.path.join(base_dir, 'datasets', 'mirbase_sequences', 'hairpin.fa.gz')
    hairpin_path_fa = os.path.join(base_dir, 'datasets', 'mirbase_sequences', 'hairpin.fa')
    hairpin_file = hairpin_path_gz if os.path.exists(hairpin_path_gz) else hairpin_path_fa
    hairpin_seqs = parse_fasta_file(hairpin_file, 'mirbase_hairpin') if os.path.exists(hairpin_file) else {}
    
    # Matching Engine across target symbols
    matched_results = {}
    matched_count = 0
    
    for rna_sym in target_names:
        key_lower = rna_sym.lower()
        seq = None
        
        # Priority 1: Mature miRNA exact match or 5p/3p arm fallback
        if key_lower in mature_seqs:
            seq = mature_seqs[key_lower]
        elif f"{key_lower}-5p" in mature_seqs:
            seq = mature_seqs[f"{key_lower}-5p"]
        elif f"{key_lower}-3p" in mature_seqs:
            seq = mature_seqs[f"{key_lower}-3p"]
        else:
            # Check if key_lower has -5p/-3p and strip it for base mature miRNA lookup
            base_mir = re.sub(r'-(5p|3p)$', '', key_lower)
            if base_mir != key_lower and base_mir in mature_seqs:
                seq = mature_seqs[base_mir]
                
        # Priority 2: GENCODE exact match (lncRNA / snoRNA / etc.)
        if not seq and key_lower in gencode_seqs:
            seq = gencode_seqs[key_lower]
            
        # Priority 3: Ensembl exact match
        if not seq and key_lower in ensembl_seqs:
            seq = ensembl_seqs[key_lower]
            
        # Priority 4: MirBase Hairpin Precursor (strictly for stem-loop precursor symbols or fallback)
        if not seq and key_lower.startswith('hsa-'):
            base_mir = re.sub(r'-(5p|3p)$', '', key_lower)
            if f"hairpin:{key_lower}" in hairpin_seqs:
                seq = hairpin_seqs[f"hairpin:{key_lower}"]
            elif f"hairpin:{base_mir}" in hairpin_seqs:
                seq = hairpin_seqs[f"hairpin:{base_mir}"]
                
        if seq and len(seq) >= 15:
            matched_results[rna_sym] = seq
            matched_count += 1
            
    print("\n" + "="*50)
    print("MULTI-KEY SEQUENCE MATCHING COMPLETE!")
    print(f"Target ncRNAs: {len(target_names)}")
    print(f"Successfully Matched Verified Sequences: {matched_count} ({matched_count/len(target_names)*100:.2f}%)")
    print("="*50)
    
    print(f"Saving matched FASTA to {out_fasta}...")
    with open(out_fasta, 'w', encoding='utf-8') as f:
        for sym, seq in matched_results.items():
            f.write(f">{sym}\n{seq}\n")
            
    print(f"Saved {len(matched_results)} verified sequence entries.")

if __name__ == "__main__":
    match_sequences()
