import click
import pandas as pd
import numpy as np
from collections import defaultdict
from ete3 import Tree
from pathlib import Path
import multiprocessing
from tqdm import tqdm
import os

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

def get_distance(tre, samples) :
    individual_pairs = defaultdict(lambda : defaultdict(dict))
    for n in tre.traverse('postorder') :
        if n.is_leaf() :
            info = n.name.split('|')
            if info[0] in samples :
                n.d = [[n.name, info[0], n.dist, 0.]]
            else :
                n.d = []
        else :
            for i, c1 in enumerate(n.children) :
                for j, c2 in enumerate(n.children[:i]) :
                    if len(c1.d) and len(c2.d) :
                        for t1 in c1.d :
                            for t2 in c2.d :
                                if (t1[1] != t2[1]) and (samples[t1[1]]['Sample_ID'] != samples[t2[1]]['Sample_ID']) :
                                    dist = t1[3] + t2[3] + (t1[2] if 'low_qual' in t1[0] else t1[2]/4) + (t2[2] if 'low_qual' in t2[0] else t2[2]/4)
                                    key = tuple(sorted([t1[1], t2[1]]))
                                    if key not in individual_pairs or individual_pairs[key] > dist :
                                        individual_pairs[key] = dist
            n.d = [[d[0], d[1], d[2], d[3] + n.dist] for c in n.children for d in c.d]
    return individual_pairs


def process_tree(data):
    nwk_file, meta_dict, threshold = data
    """Process a single tree file and return pairwise sharing results"""
    tree_results = []
    
    try:
        tree = Tree(str(nwk_file), format=1)
    except:
        return tree_results  # Skip malformed trees
    
    pairs = get_distance(tree, meta_dict)
    for (id_i, id_j), dist in pairs.items():
        shared = 1 if dist <= threshold else 0
        tree_results.append({
            'sample1': id_i,
            'sample2': id_j,
            'distance': dist,
            'shared': shared,
            'tree': nwk_file.rsplit('/', 2)[-2]  # Extract tree name from path
        })
    
    return tree_results

@click.command()
@click.option("-m",'--metadata', required=True, help='Metadata file (TSV format)')
@click.option("-t",'--tree_list', required=True, help='list of tree files (NWK format)')
@click.option('--threshold', default=0.001, type=float, 
              help='Distance threshold for strain sharing')
@click.option("-o",'--output', required=True, help='Output file for pairwise results (TSV)')
@click.option('--workers', default=8, type=int, 
              help='Number of parallel workers for processing trees')
def main(metadata, tree_list, threshold, output, workers):
    """
    Calculate strain sharing for each sample pair across multiple phylogenetic trees
    and categorize pairs using metadata.
    """
    # Load metadata
    meta_df = pd.read_csv(metadata, sep='\t').set_index('ID')
    meta_df = meta_df[~meta_df.index.duplicated(keep='first')]
    meta_dict = meta_df.to_dict(orient='index')
    
    # Get all NWK files
    with open(tree_list, 'r') as f:
        nwk_files = f.read().strip().splitlines()
    
    print(f"Found {len(nwk_files)} tree files for processing")
    
    # Setup parallel processing
    pool = multiprocessing.Pool(workers)
    tasks = [(nwk_file, meta_dict, threshold) for nwk_file in nwk_files]
    
    # Process trees in parallel
    all_results = []
    with tqdm(total=len(tasks), desc="Processing trees") as pbar:
        for tree_results in map(process_tree, tasks):
            all_results.extend(tree_results)
            pbar.update(1)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Add pair identifier (sorted)
    results_df['pair'] = results_df.apply(
        lambda row: tuple(sorted([row['sample1'], row['sample2']])), 
        axis=1
    )
    
    # Aggregate results across trees
    agg_df = results_df.groupby('pair').agg(
        sample1=('sample1', 'first'),
        sample2=('sample2', 'first'),
        trees_observed=('tree', 'count'),
        trees_shared=('shared', 'sum'),
        mean_distance=('distance', 'mean')
    ).reset_index(drop=True)
    
    # Calculate sharing rate
    agg_df['sharing_rate'] = agg_df['trees_shared'] / agg_df['trees_observed']
    
    # Add category
    agg_df['category'] = agg_df.apply(
        lambda row: categorize_pair(meta_dict[row['sample1']], meta_dict[row['sample2']]), 
        axis=1
    )
    
    # Add donor and individual information
    def get_meta_info(row, field):
        return '|'.join(sorted([meta_dict[row['sample1']].get(field, ''), meta_dict[row['sample2']].get(field, '')]))
    
    for field in ['donor', 'individual', 'Disease type']:
        agg_df[field] = agg_df.apply(lambda row: get_meta_info(row, field), axis=1)
    
    # Reorder columns
    final_df = agg_df[[
        'sample1', 'sample2', 'category', 
        'donor', 'individual', 'Disease type',
        'trees_observed', 'trees_shared', 
        'sharing_rate', 'mean_distance'
    ]]
    
    # Save results
    final_df.to_csv(output, sep='\t', index=False)
    print(f"Processed {len(nwk_files)} trees")
    print(f"Saved pairwise strain sharing results to {output}")
    
    # Print summary statistics
    print("\nSummary statistics by category:")
    summary = final_df.groupby('category').agg(
        pairs=('sharing_rate', 'count'),
        mean_sharing_rate=('sharing_rate', 'mean'),
        mean_distance=('mean_distance', 'mean')
    ).reset_index()
    print(summary.to_string(index=False))

if __name__ == '__main__':
    main()