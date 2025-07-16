import pandas as pd
import json
import click

@click.command()
@click.option('--input', '-i', required=True, help="Input TSV file from the tree analysis")
@click.option('--output', '-o', required=True, help="Output feature table TSV file")
def main(input, output):
    # Read the input TSV
    df = pd.read_csv(input, sep='\t', header=None, names=[
        'tree_path', 'cmh_pvalue', 'health_fisher', 'cohort_fisher', 
        'node', 'disease_counts', 'cohort_counts', 'ingroup', 'outgroup'
    ])

    # Filter significant rows (CMH p-value < 0.05)
    significant = df[df['cmh_pvalue'] < 0.05].copy()
    
    # Extract species name from tree_path (e.g., 'all_cohorts/Alistipes_senegalensis/uscg.nwk')
    significant['species'] = significant['tree_path'].apply(lambda x: x.split('/')[1])
    
    # Initialize a dictionary to collect sample assignments and a set for all samples
    species_genotypes = {}
    all_samples = set()

    # Process each significant species
    for _, row in significant.iterrows():
        species = row['species']
        d_counts = json.loads(row['disease_counts'])  # [[in_H, in_D], [out_H, out_D]]
        in_samples = json.loads(row['ingroup'])
        out_samples = json.loads(row['outgroup'])
        
        # Remove any sample present in both groups
        common = set(in_samples) & set(out_samples)
        if common:
            in_samples = [s for s in in_samples if s not in common]
            out_samples = [s for s in out_samples if s not in common]
        
        # Determine genotype based on disease association
        in_H, in_D = d_counts[0]
        out_H, out_D = d_counts[1]
        p_in = in_D / (in_H + in_D) if (in_H + in_D) > 0 else 0
        p_out = out_D / (out_H + out_D) if (out_H + out_D) > 0 else 0
        
        if p_in > p_out:
            geno_in = "D"
            geno_out = "H"
        else:
            geno_in = "H"
            geno_out = "D"
        
        # Store assignments for this species
        if species not in species_genotypes:
            species_genotypes[species] = {}
        for s in in_samples:
            species_genotypes[species][s] = geno_in
        for s in out_samples:
            species_genotypes[species][s] = geno_out
        
        # Track all unique samples
        all_samples.update(in_samples + out_samples)
    
    # Create feature table (species x samples)
    samples_sorted = sorted(all_samples)
    species_sorted = sorted(species_genotypes.keys())
    
    # Initialize DataFrame with '-' for missing data
    feature_df = pd.DataFrame('-', index=species_sorted, columns=samples_sorted, dtype=str)
    
    # Fill in genotypes
    for species in species_sorted:
        for sample, geno in species_genotypes[species].items():
            if sample in feature_df.columns:
                feature_df.at[species, sample] = geno
    
    # Save to TSV
    feature_df.to_csv(output, sep='\t')

if __name__ == '__main__':
    main()