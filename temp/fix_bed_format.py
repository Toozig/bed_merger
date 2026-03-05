#!/usr/bin/env python3
"""
Fix BED file formatting issues where long fields wrap to new lines.
This script combines wrapped lines back into single BED records.
"""

import sys
import re

def fix_bed_format(input_file: str, output_file: str) -> None:
    """
    Fix BED file by combining wrapped lines back into single records.
    
    Args:
        input_file: Path to input BED file with formatting issues
        output_file: Path to output fixed BED file
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_line = ""
        
        for line in infile:
            line = line.rstrip('\n\r')
            
            # If line starts with chr, it's a new BED record
            if line.startswith('chr'):
                # Write previous complete record if exists
                if current_line:
                    outfile.write(current_line + '\n')
                current_line = line
            else:
                # This is a continuation line, append to current record
                if current_line:
                    current_line += line
        
        # Write the last record
        if current_line:
            outfile.write(current_line + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_bed_format.py <input_bed> <output_bed>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"Fixing BED format: {input_file} -> {output_file}")
    fix_bed_format(input_file, output_file)
    print("Done!")
