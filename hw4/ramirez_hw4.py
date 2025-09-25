#!/home/randhal/proj/genomic_env/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def snpnexus(df, output_name):
    # Load ALL_Normal.csv
    #df = pd.read_csv("ALL_Normal.csv")
    
    # Clean chromosome names and prefix with "Chromosome"
    df["chrom"] = "Chromosome " + df["chrom"].astype(str).str.replace("chr", "", case=False)
    
    # Ensure positions are integers
    df["left"] = df["left"].astype(int)
    
    # Add a sample ID column (all "1")
    df["sample"] = 1
    
    # Arrange columns: Chromosome, Position, Ref, Alt, Sample
    snpnexus_df = df[["chrom", "left", "ref_seq", "alt_seq", "sample"]]
    
    # Save as tab-delimited
    snpnexus_df.to_csv(output_name, sep="\t", index=False, header=False)

if __name__=='__main__':
    data_normal = pd.read_csv("ALL_Normal.csv")
    data_tumor = pd.read_csv("ALL_Tumor.csv")
    
    snpnexus(data_normal, "data_normal.txt")
    snpnexus(data_tumor, "data_tumor.txt")
