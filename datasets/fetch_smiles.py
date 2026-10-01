import pandas as pd
import requests
import time
import os
import glob

excel_files = glob.glob('d:/fydp/datasets/ncrnadrug/*.xlsx')
all_drugs = set()

for f in excel_files:
    try:
        df = pd.read_excel(f)
        drug_col = None
        for col in df.columns:
            if 'drug' in col.lower() or 'compound' in col.lower():
                drug_col = col
                break
        
        if drug_col:
            drugs = df[drug_col].dropna().unique()
            all_drugs.update(drugs)
    except:
        pass

results = []
for i, drug in enumerate(all_drugs):
    drug_str = str(drug).strip()
    if not drug_str: continue
        
    try:
        # Request both SMILES and IsomericSMILES just in case
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_str}/property/IsomericSMILES,CanonicalSMILES/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                props = data['PropertyTable']['Properties'][0]
                smiles = props.get('IsomericSMILES', props.get('CanonicalSMILES', props.get('SMILES', '')))
                if smiles:
                    results.append({'Drug_Name': drug_str, 'SMILES': smiles, 'Found': True})
                else:
                    results.append({'Drug_Name': drug_str, 'SMILES': '', 'Found': False})
            else:
                results.append({'Drug_Name': drug_str, 'SMILES': '', 'Found': False})
        else:
            results.append({'Drug_Name': drug_str, 'SMILES': '', 'Found': False})
            
    except:
        results.append({'Drug_Name': drug_str, 'SMILES': '', 'Found': False})
    
    time.sleep(0.2)

out_file = 'd:/fydp/datasets/drugbank/drug_smiles.csv'
out_df = pd.DataFrame(results)
out_df.to_csv(out_file, index=False)
print(f"Fixed and saved {len(out_df[out_df['Found'] == True])} valid SMILES to {out_file}")
