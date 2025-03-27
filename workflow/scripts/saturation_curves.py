#!/usr/bin/env python3

import os
import subprocess
import random
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate sequencing depth saturation curves for nanoCUT&Tag data.')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Root directory containing all samples')
    parser.add_argument('--samples', type=str, nargs='+', 
                        help='Sample IDs to process. If not provided, all samples in the output_dir will be processed.')
    parser.add_argument('--modalities', type=str, nargs='+', default=['H3K27ac', 'H3K27me3'],
                        help='Modalities to analyze for each sample (default: H3K27ac and H3K27me3)')
    parser.add_argument('--downsample_ratios', type=float, nargs='+', 
                        default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                        help='Ratios for downsampling BAM files')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads to use for parallel processing')
    parser.add_argument('--output_file', type=str, default='saturation_curves.pdf',
                        help='Output file for the saturation curves plot')
    parser.add_argument('--temp_dir', type=str, default='temp',
                        help='Directory for temporary files')
    return parser.parse_args()

def get_all_samples(output_dir):
    """Get all sample directories in the output directory"""
    return [d for d in os.listdir(output_dir) 
            if os.path.isdir(os.path.join(output_dir, d)) and not d.startswith('.')]

def downsample_bam(bam_file, output_bam, ratio, seed=42):
    """Downsample a BAM file using samtools"""
    print(f"Downsampling {bam_file} to {ratio:.1%}")
    
    # Convert ratio to proper samtools format (needs to be formatted as SEED.FRACTION)
    # For samtools, the fraction needs to be an integer between 0 and 1, where 0.1 = 10%
    # The correct format is: "seed + decimal point + integer fraction"
    fraction = str(ratio).split('.')[-1]
    if fraction == '0':
        fraction = '1'  # Handle 1.0 -> use fraction '1'
    
    cmd = f"samtools view -b -s {seed}.{fraction} {bam_file} > {output_bam}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Index the downsampled BAM
    cmd = f"samtools index {output_bam}"
    subprocess.run(cmd, shell=True, check=True)
    return output_bam

def generate_fragments_from_bam(bam_file, output_fragments):
    """Generate fragments file from BAM file using sinto and ensure it's properly sorted"""
    print(f"Generating fragments from {bam_file}")
    # Get directory for output fragments
    output_dir = os.path.dirname(output_fragments)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Run sinto to generate fragments
    temp_fragments = f"{output_fragments}.tmp"
    cmd = f"sinto fragments -b {bam_file} -f {temp_fragments} --barcode_regex '[^:]*'"
    subprocess.run(cmd, shell=True, check=True)
    
    # Sort the fragments file (important for tabix indexing)
    print(f"Sorting fragments file")
    # Sort by chromosome and then by start position
    sorted_fragments = output_fragments
    
    # Use awk to handle potential chromosome naming issues and sort numerically
    sort_cmd = f"""awk '{{print $0}}' {temp_fragments} | sort -k1,1 -k2,2n > {sorted_fragments}"""
    subprocess.run(sort_cmd, shell=True, check=True)
    
    # Remove temporary file
    os.remove(temp_fragments)
    
    # Compress the sorted fragments file
    cmd = f"bgzip -f {sorted_fragments}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Index the compressed fragments file
    cmd = f"tabix -p bed {sorted_fragments}.gz"
    subprocess.run(cmd, shell=True, check=True)
    
    return f"{sorted_fragments}.gz"

def calculate_metrics(fragments_file):
    """Calculate median unique fragments per cell and mean reads per cell from fragments file"""
    print(f"Calculating metrics for {fragments_file}")
    
    # Read the fragments file (it's a TSV with columns: chr, start, end, barcode, count)
    try:
        df = pd.read_csv(fragments_file, sep='\t', header=None, 
                        names=['chr', 'start', 'end', 'barcode', 'count'],
                        compression='gzip')
    except Exception as e:
        print(f"Error reading fragments file: {e}")
        # Try to read without assuming compression
        try:
            if fragments_file.endswith('.gz'):
                # Try reading the file directly with pandas' built-in gzip support
                df = pd.read_csv(fragments_file, sep='\t', header=None, 
                            names=['chr', 'start', 'end', 'barcode', 'count'],
                            compression='gzip')
            else:
                # Try reading uncompressed
                df = pd.read_csv(fragments_file, sep='\t', header=None, 
                            names=['chr', 'start', 'end', 'barcode', 'count'])
        except Exception as e2:
            print(f"Still can't read file after trying alternatives: {e2}")
            return None
    
    # Group by barcode to calculate fragments per cell
    fragments_per_cell = df.groupby('barcode').size()
    
    # Calculate median unique fragments per cell
    median_fragments = fragments_per_cell.median()
    
    # Calculate total reads across all cells (sum of counts)
    total_reads = df['count'].sum()
    
    # Calculate mean reads per cell
    num_cells = len(fragments_per_cell)
    mean_reads_per_cell = total_reads / num_cells if num_cells > 0 else 0
    
    return {
        'median_fragments_per_cell': median_fragments,
        'mean_reads_per_cell': mean_reads_per_cell,
        'num_cells': num_cells
    }

def process_sample_modality(args, sample, modality, downsample_ratio):
    """Process a single sample and modality at a specific downsample ratio"""
    try:
        sample_dir = os.path.join(args.output_dir, sample, modality)
        bam_file = os.path.join(sample_dir, 'possorted_bam.bam')
        
        if not os.path.exists(bam_file):
            print(f"BAM file not found for {sample}/{modality}, skipping")
            return None
        
        # Create temporary directory for downsampled files
        temp_dir = os.path.join(args.temp_dir, sample, modality)
        os.makedirs(temp_dir, exist_ok=True)
        
        # Define output files
        downsampled_bam = os.path.join(temp_dir, f"downsampled_{downsample_ratio:.1f}.bam")
        fragments_file = os.path.join(temp_dir, f"fragments_{downsample_ratio:.1f}.tsv")
        
        # Downsample BAM file
        if downsample_ratio < 1.0:
            downsample_bam(bam_file, downsampled_bam, downsample_ratio)
            bam_to_use = downsampled_bam
        else:
            bam_to_use = bam_file
        
        # Generate fragments from BAM
        fragments_path = generate_fragments_from_bam(bam_to_use, fragments_file)
        
        # Calculate metrics
        metrics = calculate_metrics(fragments_path)
        
        if metrics:
            # Add sample information
            metrics['sample'] = sample
            metrics['modality'] = modality
            metrics['downsample_ratio'] = downsample_ratio
            
            try:
                if downsample_ratio < 1.0:
                    os.remove(downsampled_bam)
                    os.remove(f"{downsampled_bam}.bai")  # Remove BAM index
                # Remove fragments files
                os.remove(fragments_path)
                os.remove(f"{fragments_path}.tbi")
            except Exception as e:
                print(f"Warning: Error deleting temporary files for {sample}/{modality} ratio {downsample_ratio}: {e}")
            
            return metrics
        return None
    except Exception as e:
        print(f"Error processing {sample}/{modality} at {downsample_ratio:.1%}: {e}")
        return None
    

def plot_saturation_curves(results_df, output_file):
    """Plot saturation curves based on the collected metrics"""
    plt.figure(figsize=(12, 8))
    
    # Define line styles for modalities
    modality_styles = {"H3K27ac": "-", "H3K27me3": "--"}
    
    # Create a color palette for different samples
    samples = results_df['sample'].unique()
    # Use a colorblind-friendly palette with distinct colors
    sample_colors = plt.cm.tab10(range(len(samples)))
    sample_palette = {sample: sample_colors[i] for i, sample in enumerate(samples)}
    
    # Plot each sample and modality with different styles
    for sample in samples:
        sample_data = results_df[results_df['sample'] == sample]
        sample_color = sample_palette[sample]
        
        for modality in sample_data['modality'].unique():
            subset = sample_data[sample_data['modality'] == modality]
            subset = subset.sort_values('downsample_ratio')  # Sort by downsampling ratio
            
            linestyle = modality_styles.get(modality, "-")
            
            # Plot with label only the first time the sample-modality combination appears
            label = f"{sample} - {modality}" if 'label' not in locals() or label != f"{sample} - {modality}" else None
            
            plt.plot(subset['mean_reads_per_cell'], subset['median_fragments_per_cell'],
                    marker='o', linestyle=linestyle, color=sample_color, alpha=0.8,
                    label=label)
    
    plt.title('Sequencing Depth Saturation Curves', fontsize=16)
    plt.xlabel('Mean Reads per Cell', fontsize=14)
    plt.ylabel('Median Unique Fragments per Cell', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Sample - Modality')
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saturation curve plot saved to {output_file}")
    
    csv_file = output_file.replace('.pdf', '.csv')
    results_df.to_csv(csv_file, index=False)
    print(f"Results data saved to {csv_file}")

def main():
    args = parse_arguments()
    
    # Create temporary directory
    os.makedirs(args.temp_dir, exist_ok=True)
    
    # Get samples to process
    samples = args.samples if args.samples else get_all_samples(args.output_dir)
    print(f"Processing samples: {', '.join(samples)}")
    
    # Prepare all tasks
    tasks = []
    for sample in samples:
        for modality in args.modalities:
            for ratio in args.downsample_ratios:
                # Check if sample/modality exists
                if os.path.exists(os.path.join(args.output_dir, sample, modality)):
                    tasks.append((sample, modality, ratio))
    
    # Process tasks sequentially for debugging or in parallel for speed
    results = []
    
    # For debugging, process sequentially
    # for sample, modality, ratio in tasks:
    #     result = process_sample_modality(args, sample, modality, ratio)
    #     if result:
    #         results.append(result)
    
    # For production, process in parallel
    with Pool(args.threads) as pool:
        process_func = partial(process_sample_modality, args)
        pool_results = pool.starmap(process_func, [(sample, modality, ratio) for sample, modality, ratio in tasks])
        results = [r for r in pool_results if r is not None]
    
    # Convert results to DataFrame
    if results:
        results_df = pd.DataFrame(results)
        
        # Plot the results
        plot_saturation_curves(results_df, args.output_file)
    else:
        print("No valid results were generated. Please check your input data and parameters.")
    
    # Clean up temporary files if needed
    import shutil
    shutil.rmtree(args.temp_dir)

if __name__ == "__main__":
    main()
