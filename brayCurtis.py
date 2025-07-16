import os
import numpy as np
import pandas as pd
import click
from itertools import combinations
from scipy.spatial.distance import braycurtis

def parse_profile(profile_path):
    """Parse microbiome profile file into DataFrame"""
    with open(profile_path, 'r') as f:
        # Parse header to extract sample IDs
        header = next(f).split()
        sample_paths = header[1:-1]  # Skip "#Species(RPKM)" and "#Taxonomy"
        sample_ids = [os.path.basename(os.path.dirname(p)) for p in sample_paths]
        
        # Read data into DataFrame
        data = []
        species_names = []
        # taxonomies = []
        
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            
            # Extract species name, abundances, and taxonomy
            taxonomy = parts[-1]
            if taxonomy.find('d__Viruses') < 0 :
                species_names.append(parts[0])
                data.append([float(x) for x in parts[2:2+len(sample_ids)]])
    
    # Create DataFrame with samples as columns
    df = pd.DataFrame(data, index=species_names, columns=sample_ids).T
    return df

def categorize_pair(row1, row2):
    """Categorize sample pair based on metadata"""
    type1 = row1['Disease type']
    type2 = row2['Disease type']
    
    if type1 == 'healthy' and type2 == 'healthy':
        if row1['individual'] == row2['individual']:
            return 'healthy_same_individual'
        return 'healthy_different_individuals'
    
    if type1 != 'healthy' and type2 != 'healthy':
        if type1 != 'before_FMT' and type2 != 'before_FMT':
            if row1['individual'] == row2['individual']:
                return 'same_patients'
            if row1['donor'] == row2['donor']:
                return 'different_patients_same_donor'
            if row1['donor'] != row2['donor']:
                return 'different_patients'
        elif row1['individual'] == row2['individual'] and (type1 != 'before_FMT' or type2 != 'before_FMT'):
            if row1['Disease type'] == 'responder' or row2['Disease type'] == 'responder':
                return 'before_after_FMT_same_responder'
            elif row1['Disease type'] == 'non-responder' or row2['Disease type'] == 'non-responder':
                return 'before_after_FMT_same_non_responder'
            # return 'before_after_FMT_same_individual'
    
    # One healthy, one patient
    if type1 != 'healthy':
        patient, healthy = row1, row2
    else:
        patient, healthy = row2, row1
    
    if patient['donor'] == healthy['individual'] and patient['Disease type'] != 'before_FMT':
        if patient['Disease type'] == 'responder':
            return 'patient_vs_its_donor_responder'
        elif patient['Disease type'] == 'non-responder':
            return 'patient_vs_its_donor_non_responder'
    
    return 'other'  # Default category for unmatched pairs

@click.command()
@click.option('-p', '--profile', required=True, help='Microbiome profile file path')
@click.option('-m', '--metadata', required=True, help='Sample metadata file path')
@click.option('-o', '--output', default='distances.tsv', help='Output file path')
def main(profile, metadata, output):
    """Calculate pairwise Bray-Curtis distances with sample categorization"""
    # Load data
    profile_df = parse_profile(profile)
    meta_df = pd.read_csv(metadata, sep='\t')
    
    # Set index for metadata
    meta_df = meta_df.set_index('ID')
    
    # Filter to samples present in both files
    common_samples = list(set(profile_df.index) & set(meta_df.index))
    profile_df = profile_df.loc[common_samples]
    meta_df = meta_df.loc[common_samples]
    meta_df = meta_df[~meta_df.index.duplicated(keep='first')]
    
    # Calculate Bray-Curtis distances
    results = []
    sample_pairs = list(combinations(profile_df.index, 2))
    
    for (s1, s2) in sample_pairs:
        dist = braycurtis(profile_df.loc[s1], profile_df.loc[s2])
        category = categorize_pair(meta_df.loc[s1], meta_df.loc[s2])
        results.append({
            'Sample1': s1,
            'Sample2': s2,
            'BrayCurtis': dist,
            'Category': category
        })
    
    # Create and save results
    result_df = pd.DataFrame(results)
    result_df.to_csv(output, sep='\t', index=False)
    print(f"Saved {len(result_df)} pairwise distances to {output}")

if __name__ == '__main__':
    main()