import pandas as pd
from scipy.stats import hypergeom
import sys

def hypergeometric_test(input_file):
    # Read the input file
    df = pd.read_csv(input_file, sep='\t')
    
    # Define the threshold for high Fst
    fst_threshold = 0.25
    
    # Prepare a list to store results
    results = []

    # Group by scaffold
    grouped = df.groupby('CHROM')
    
    # Total number of SNPs in the dataset
    total_snps = len(df)
    
    for scaffold, group in grouped:
        # Count total SNPs in the scaffold
        total_snps_scaffold = len(group)
        
        # Count high Fst SNPs in the scaffold
        high_fst_snps = group[group['WEIR_AND_COCKERHAM_FST'] > fst_threshold]
        num_high_fst_scaffold = len(high_fst_snps)
        
        # Count total high Fst SNPs in the entire dataset
        total_high_fst = len(df[df['WEIR_AND_COCKERHAM_FST'] > fst_threshold])
        
        # Perform hypergeometric test
        # k = number of high Fst SNPs in the scaffold
        # M = total number of high Fst SNPs in the population
        # n = total number of SNPs in the scaffold
        # N = total number of SNPs in the population
        p_value = hypergeom.sf(num_high_fst_scaffold - 1, total_snps, total_high_fst, total_snps_scaffold)
        
        # Append the result
        results.append((scaffold, p_value, total_snps_scaffold, num_high_fst_scaffold))
    
    # Create a DataFrame for results
    results_df = pd.DataFrame(results, columns=['Scaffold', 'p-value', 'total_snps_scaffold', 'num_high_fst_snps_scaffold'])
    
    # Output the results to a tab-delimited file
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f'Results saved to {output_file}')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 HypergeometricTest.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    hypergeometric_test(input_file)
