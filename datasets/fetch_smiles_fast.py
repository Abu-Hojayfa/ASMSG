import pandas as pd
import requests
import time
import os

old_df = pd.read_csv('d:/fydp/datasets/drugbank/drug_smiles.csv')
drugs_to_fetch = old_df[old_df['Found'] == True]['Drug_Name'].tolist()

print(f"Fetching SMILES for {len(drugs_to_fetch)} drugs...")
results = []

for i, drug in enumerate(drugs_to_fetch):
    drug_str = str(drug).strip()
    if i % 100 == 0:
        print(f"[{i}/{len(drugs_to_fetch)}] Processing {drug_str}...")
        
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_str}/property/IsomericSMILES,CanonicalSMILES/JSON"
        response = requests.get(url, timeout=5)
        
        smiles = ''
        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                props = data['PropertyTable']['Properties'][0]
                smiles = props.get('IsomericSMILES', props.get('CanonicalSMILES', props.get('SMILES', '')))
        
        results.append({'Drug_Name': drug_str, 'SMILES': smiles, 'Found': True if smiles else False})
            
    except Exception as e:
        results.append({'Drug_Name': drug_str, 'SMILES': '', 'Found': False})
    
    time.sleep(0.2)

# Add back the ones that were False originally
false_drugs = old_df[old_df['Found'] == False]['Drug_Name'].tolist()
for drug in false_drugs:
    results.append({'Drug_Name': str(drug).strip(), 'SMILES': '', 'Found': False})

out_file = 'd:/fydp/datasets/drugbank/drug_smiles.csv'
out_df = pd.DataFrame(results)
out_df.to_csv(out_file, index=False)
print(f"\nDone! Overwrote {out_file} with actual SMILES strings.")
