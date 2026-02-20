import pandas as pd
from pathlib import Path

# Read all stratification files
dfs = []
for tsv in snakemake.input:
    stratification = Path(tsv).stem
    df = pd.read_csv(
        tsv, 
        sep="\t", 
        names=["depth", stratification], 
        dtype={"depth": int, stratification: int})
    df.set_index("depth", inplace=True)
    dfs.append(df)

# Concatenate all stratifications into one dataframe
df = pd.concat(dfs, axis=1)
df = df.fillna(0)

# Save the final dataframe
df.to_csv(snakemake.output.tsv, sep="\t")

# Normalize the data
df = df.div(df.sum(axis=0), axis=1)
df.to_csv(snakemake.output.tsv2, sep="\t")
