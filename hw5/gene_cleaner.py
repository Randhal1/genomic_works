#!/Users/randhal/Documents/Interpreters/xptopics/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Drop the duplicates of the AML genes file
def clean_genes(data):
    df = pd.read_csv(data)
    aa = df["SYMBOL"].drop_duplicates()
    aa.to_csv("AML_gene_clean.csv", index=False)


if __name__=='__main__':
    clean_genes("AML_gene.csv")
