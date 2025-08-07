#!/usr/bin/env python3
import os
import glob
import shutil

def main():
    # Working directory (adjust if needed)
    directory = "."
    
    # Source file to copy from
    src_file = os.path.join(directory, "Winter24Run3_V1_MC_L2Relative_AK4PUPPI.txt")
    if not os.path.isfile(src_file):
        print(f"Source file not found: {src_file}")
        return

    # Pattern for files that indicate an L2L3Residual correction for DATA
    pattern = os.path.join(directory, "Prompt24_*_DATA_L2L3Residual_AK4PFPuppi.txt")
    matching_files = glob.glob(pattern)
    
    if not matching_files:
        print(f"No matching files found with pattern: {pattern}")
        return

    for old_file in matching_files:
        # Construct new file name by replacing "L2L3Residual" with "L2Relative"
        base_name = os.path.basename(old_file)
        new_base = base_name.replace("L2L3Residual", "L2Relative")
        new_file = os.path.join(directory, new_base)

        # Check if destination file already exists.
        if os.path.exists(new_file):
            print(f"Destination file {new_file} already exists. Skipping copy.")
            continue
        
        # Copy the source file to the new destination
        shutil.copy2(src_file, new_file)
        print(f"Copied {src_file} to {new_file}")

if __name__ == '__main__':
    main()
