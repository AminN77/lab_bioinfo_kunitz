import pandas as pd

df = pd.read_csv("pdb_kunitz_raw.csv", skiprows=1)
df.columns = ["entity_id","pdb_id","resolution","sequence","chain","entry_id","extra"]
df[["entity_id","pdb_id","resolution","sequence","entry_id"]] = \
    df[["entity_id","pdb_id","resolution","sequence","entry_id"]].ffill()
df = df.dropna(subset=["chain"]).drop_duplicates(subset=["pdb_id","chain"])

with open("kunitz_pdb.fasta","w") as f:
    for _, r in df.iterrows():
        f.write(f">{r.pdb_id}_{r.chain}\n{r.sequence}\n")

print(f"Chains kept: {len(df)}")
