#!/usr/bin/env python3
"""
BED File Analysis and Excel Report Generator

This script generates comprehensive Excel reports for BED file analysis including:
- Individual sheets for each BED file
- Statistics summary with coding region overlap analysis
- GTF-based coding region extraction
"""

import pandas as pd
import sys
import os
import subprocess
import tempfile
from pathlib import Path
import argparse

try:
    import gtfparse
except ImportError:
    print("Warning: gtfparse not available. Will use basic GTF parsing.")
    gtfparse = None


def extract_coding_regions_from_gtf(gtf_file, output_bed):
    """Extract coding regions from GTF file and save as BED format"""
    print(f"Extracting coding regions from {gtf_file}...")
    
    if gtfparse:
        # Use gtfparse for better GTF handling
        try:
            df = gtfparse.read_gtf(gtf_file)
            
            # Filter for CDS features
            cds_df = df[df['feature'] == 'CDS'].copy()
            
            if len(cds_df) == 0:
                print("No CDS features found in GTF file")
                return 0, 0
            
            # Convert to BED format (0-based coordinates)
            bed_df = pd.DataFrame({
                'chr': cds_df['seqname'],
                'start': cds_df['start'] - 1,  # Convert to 0-based
                'end': cds_df['end']
            })
            
            # Sort and merge overlapping regions
            bed_df = bed_df.sort_values(['chr', 'start'])
            bed_df.to_csv(output_bed + '.tmp', sep='\t', header=False, index=False)
            
            # Use bedtools to merge overlapping regions
            subprocess.run([
                'bedtools', 'merge', '-i', output_bed + '.tmp'
            ], stdout=open(output_bed, 'w'), check=True)
            
            # Clean up temp file
            os.unlink(output_bed + '.tmp')
            
        except Exception as e:
            print(f"Error using gtfparse: {e}. Falling back to basic parsing.")
            # Use basic parsing instead
            use_basic_parsing = True
        else:
            use_basic_parsing = False
    else:
        use_basic_parsing = True
    
    if use_basic_parsing:
        # Fallback to basic awk-based extraction
        subprocess.run([
            'bash', '-c',
            f"awk -F'\\t' '$3 == \"CDS\" {{print $1 \"\\t\" ($4-1) \"\\t\" $5}}' {gtf_file} | "
            f"sort -k1,1 -k2,2n | bedtools merge -i - > {output_bed}"
        ], check=True)
    
    # Calculate statistics
    coding_count = int(subprocess.run(['wc', '-l', output_bed], 
                                    capture_output=True, text=True).stdout.split()[0])
    
    coding_bp = int(subprocess.run([
        'awk', '-F\\t', '{sum += $3 - $2} END {print sum}', output_bed
    ], capture_output=True, text=True).stdout.strip())
    
    print(f"Extracted {coding_count} coding regions ({coding_bp:,} bp total)")
    return coding_count, coding_bp


def calculate_bed_stats(bed_file):
    """Calculate comprehensive statistics for a BED file"""
    try:
        # Read BED file with more robust handling
        with open(bed_file, 'r', newline='') as f:
            # Read and clean lines to handle different line endings
            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        # Parse lines into DataFrame
        data = []
        for line in lines:
            # Try tab first, then space separation
            parts = line.split('\t')
            if len(parts) < 3:
                parts = line.split()
            
            if len(parts) >= 3:  # Valid BED line needs at least 3 columns
                data.append(parts)
        
        if not data:
            print(f"No valid BED data found in {bed_file}")
            return None, None
        
        df = pd.DataFrame(data)
        
        # Handle different numbers of columns
        column_names = ['chr', 'start', 'end']
        if df.shape[1] > 3:
            column_names.extend([f'col_{i}' for i in range(4, df.shape[1] + 1)])
        
        df.columns = column_names[:df.shape[1]]
        
        # Coerce coordinate columns to numeric and drop invalid rows
        if 'start' in df.columns and 'end' in df.columns:
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
            df['end'] = pd.to_numeric(df['end'], errors='coerce')
            df = df.dropna(subset=['start', 'end'])
            # Ensure integer dtype for calculations and writing
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            # Filter out malformed intervals
            df = df[df['end'] > df['start']]
        
        # Calculate segment lengths
        lengths = df['end'] - df['start']
        total_bp = lengths.sum()
        
        stats = {
            'n_segments': len(df),
            'total_bp': int(total_bp),
            'min_length': int(lengths.min()),
            'max_length': int(lengths.max()),
            'mean_length': float(lengths.mean()),
            'median_length': float(lengths.median())
        }
        
        return stats, df
        
    except Exception as e:
        print(f"Error processing {bed_file}: {e}")
        return None, None


def calculate_coding_overlap(bed_df, coding_bed_file):
    """Calculate overlap between BED regions and coding regions"""
    try:
        # Create temporary BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
            bed_df[['chr', 'start', 'end']].to_csv(tmp, sep='\t', header=False, index=False)
            tmp_bed = tmp.name
        
        # Count segments that overlap with coding regions
        result = subprocess.run([
            'bedtools', 'intersect', '-a', tmp_bed, '-b', coding_bed_file, '-u'
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            intersecting_count = len([line for line in result.stdout.strip().split('\n') 
                                    if line.strip()]) if result.stdout.strip() else 0
        else:
            intersecting_count = 0
        
        # Calculate base pair overlap
        result_bp = subprocess.run([
            'bedtools', 'intersect', '-a', tmp_bed, '-b', coding_bed_file
        ], capture_output=True, text=True)
        
        if result_bp.returncode == 0 and result_bp.stdout.strip():
            overlap_lines = [line for line in result_bp.stdout.strip().split('\n') if line.strip()]
            overlap_bp = sum(int(line.split('\t')[2]) - int(line.split('\t')[1]) 
                           for line in overlap_lines)
        else:
            overlap_bp = 0
        
        # Clean up
        os.unlink(tmp_bed)
        
        return intersecting_count, overlap_bp
        
    except Exception as e:
        print(f"Error calculating coding overlap: {e}")
        return 0, 0


def generate_excel_report(bed_files, coding_bed_file, output_excel):
    """Generate comprehensive Excel report with multiple sheets"""
    print(f"Generating Excel report: {output_excel}")
    
    all_stats = []
    
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        
        for bed_file in bed_files:
            file_name = os.path.basename(bed_file)
            print(f"Processing {file_name}...")
            
            # Calculate basic statistics
            stats, bed_df = calculate_bed_stats(bed_file)
            if stats is None:
                print(f"Skipping {file_name} due to errors")
                continue
            
            # Calculate coding overlap
            coding_segments, coding_bp = calculate_coding_overlap(bed_df, coding_bed_file)
            coding_percent = (coding_bp / stats['total_bp'] * 100) if stats['total_bp'] > 0 else 0
            
            # Add additional stats
            stats.update({
                'file_name': file_name,
                'coding_segments': coding_segments,
                'coding_bp': coding_bp,
                'coding_percent': coding_percent
            })
            
            all_stats.append(stats)
            
            # Create sheet for this BED file
            sheet_name = file_name.replace('.bed', '')[:31]  # Excel sheet name limit
            
            # Add column headers if missing
            if 'chr' in bed_df.columns:
                bed_df.to_excel(writer, sheet_name=sheet_name, index=False)
            else:
                # Add standard BED headers
                bed_df_copy = bed_df.copy()
                bed_df_copy.columns = ['chr', 'start', 'end'] + [f'col_{i}' for i in range(4, len(bed_df_copy.columns) + 1)]
                bed_df_copy.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Create statistics summary sheet
        if all_stats:
            stats_df = pd.DataFrame(all_stats)
            
            # Reorder columns for better readability
            column_order = ['file_name', 'n_segments', 'total_bp', 'min_length', 'max_length', 
                          'mean_length', 'median_length', 'coding_segments', 'coding_bp', 'coding_percent']
            
            stats_df = stats_df[column_order]
            
            # Round numeric columns
            numeric_cols = ['mean_length', 'median_length', 'coding_percent']
            for col in numeric_cols:
                if col in stats_df.columns:
                    stats_df[col] = stats_df[col].round(2)
            
            # Format column names for better readability
            stats_df.columns = ['File Name', 'Segments', 'Total BP', 'Min Length', 'Max Length',
                              'Mean Length', 'Median Length', 'Coding Segments', 'Coding BP', 'Coding %']
            
            stats_df.to_excel(writer, sheet_name='Statistics_Summary', index=False)
            
            print(f"\nExcel report generated successfully!")
            print(f"Processed {len(all_stats)} BED files")
            print(f"Report saved to: {output_excel}")
            
            # Print summary statistics
            print(f"\nSummary:")
            print(f"  Total segments across all files: {stats_df['Segments'].sum():,}")
            print(f"  Total base pairs: {stats_df['Total BP'].sum():,}")
            print(f"  Average coding overlap: {stats_df['Coding %'].mean():.2f}%")
        else:
            print("No valid BED files processed!")


def main():
    parser = argparse.ArgumentParser(description='Generate Excel report for BED file analysis')
    parser.add_argument('bed_files', nargs='+', help='BED files to analyze')
    parser.add_argument('--gtf', required=True, help='GTF file for coding region extraction')
    parser.add_argument('--output', required=True, help='Output Excel file')
    parser.add_argument('--coding-bed', help='Pre-generated coding regions BED file (optional)')
    
    args = parser.parse_args()
    
    # Check input files
    for bed_file in args.bed_files:
        if not os.path.exists(bed_file):
            print(f"Error: BED file not found: {bed_file}")
            sys.exit(1)
    
    if not os.path.exists(args.gtf):
        print(f"Error: GTF file not found: {args.gtf}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir:  # Only create if there's actually a directory path
        os.makedirs(output_dir, exist_ok=True)
    
    # Extract coding regions if not provided
    if args.coding_bed:
        coding_bed_file = args.coding_bed
        # Ensure directory exists for coding bed file
        coding_bed_dir = os.path.dirname(coding_bed_file)
        if coding_bed_dir:
            os.makedirs(coding_bed_dir, exist_ok=True)
        # If coding bed file doesn't exist, create it
        if not os.path.exists(coding_bed_file):
            extract_coding_regions_from_gtf(args.gtf, coding_bed_file)
    else:
        coding_bed_file = args.output.replace('.xlsx', '_coding_regions.bed')
        extract_coding_regions_from_gtf(args.gtf, coding_bed_file)
    
    # Generate Excel report
    generate_excel_report(args.bed_files, coding_bed_file, args.output)


if __name__ == "__main__":
    main()
