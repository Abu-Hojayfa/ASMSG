import pandas as pd
import os

def explore_rnadisease():
    file_path = r"d:\fydp\datasets\rnadisease_v4\alldata.xlsx"
    print(f"Loading {file_path} (This might take 10-20 seconds because it's a large Excel file)...")
    
    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"Error loading file: {e}")
        return
    
    print("\n" + "="*50)
    print("Columns found in RNADisease v4.0:")
    print("="*50)
    for col in df.columns:
        print(f" - {col}")
        
    print("\n" + "="*50)
    print("First 3 rows of data:")
    print("="*50)
    print(df.head(3).to_string())

if __name__ == "__main__":
    explore_rnadisease()
