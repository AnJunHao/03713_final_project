from pathlib import Path
import subprocess
from typing import NamedTuple
from pipeline.monitor import monitor_jobs
from pipeline.bedtool_preprocess import BedtoolConfig, load_bedtool_config
from tabulate import tabulate

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {species} {tissue} conserved regions analysis"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Process conserved regions
echo "[STEP 1] Processing conserved regions"
if [ "{direction}" = "{species_1}_to_{species_2}" ]; then
    # Process species_1 to species_2 conserved regions
    bedtools closest -a {conserved_bed} \\
                     -b {tss_file} \\
                     -d > {output_dir}/{direction}_{tissue}_conserved.TSS.bed

    # Step 2: Classify as promoters vs enhancers
    echo "[STEP 2] Classifying conserved regions as enhancers and promoters"
    awk '$NF <= 5000' {output_dir}/{direction}_{tissue}_conserved.TSS.bed > {output_dir}/{direction}_{tissue}_conserved_promoters.bed
    awk '$NF > 5000' {output_dir}/{direction}_{tissue}_conserved.TSS.bed > {output_dir}/{direction}_{tissue}_conserved_enhancers.bed
elif [ "{direction}" = "{species_2}_to_{species_1}" ]; then
    # Process species_2 to species_1 conserved regions
    bedtools closest -a {conserved_bed} \\
                     -b {tss_file} \\
                     -d > {output_dir}/{direction}_{tissue}_conserved.TSS.bed

    # Step 2: Classify as promoters vs enhancers
    echo "[STEP 2] Classifying conserved regions as enhancers and promoters"
    awk '$NF <= 5000' {output_dir}/{direction}_{tissue}_conserved.TSS.bed > {output_dir}/{direction}_{tissue}_conserved_promoters.bed
    awk '$NF > 5000' {output_dir}/{direction}_{tissue}_conserved.TSS.bed > {output_dir}/{direction}_{tissue}_conserved_enhancers.bed
fi

# Step 3: Count peaks
echo "[STEP 3] Counting conserved regions"
total_regions=$(wc -l {conserved_bed} | awk '{{print $1}}')
promoters=$(wc -l {output_dir}/{direction}_{tissue}_conserved_promoters.bed | awk '{{print $1}}')
enhancers=$(wc -l {output_dir}/{direction}_{tissue}_conserved_enhancers.bed | awk '{{print $1}}')

echo "Total {direction} {tissue} conserved regions: $total_regions"
echo "Conserved promoters: $promoters"
echo "Conserved enhancers: $enhancers"

echo "Job finished"
"""

shared_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing shared across-species analysis for {tissue}"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# First we need to ensure we have the target species promoters and enhancers
echo "[STEP 0] Preparing target species {species_2} {tissue} promoters and enhancers"
bedtools closest -a {species_2_peak_file} \\
                 -b {species_2_tss_file} \\
                 -d > {output_dir}/{species_2}_{tissue}_TSS_annotated.bed

awk '$NF <= 5000' {output_dir}/{species_2}_{tissue}_TSS_annotated.bed > {output_dir}/{species_2}_{tissue}_promoters.bed
awk '$NF > 5000'  {output_dir}/{species_2}_{tissue}_TSS_annotated.bed > {output_dir}/{species_2}_{tissue}_enhancers.bed

# Step 1: Find shared promoters across species
echo "[STEP 1] Finding shared promoters across species"
bedtools intersect -a {species_1_to_species_2_promoters} \\
                  -b {output_dir}/{species_2}_{tissue}_promoters.bed \\
                  -u > {output_dir}/{tissue}_promoters_shared_across_species.bed

# Step 2: Find shared enhancers across species 
echo "[STEP 2] Finding shared enhancers across species"
bedtools intersect -a {species_1_to_species_2_enhancers} \\
                  -b {output_dir}/{species_2}_{tissue}_enhancers.bed \\
                  -u > {output_dir}/{tissue}_enhancers_shared_across_species.bed

# Step 3: Count shared elements
echo "[STEP 3] Counting shared elements across species"
shared_promoters=$(wc -l {output_dir}/{tissue}_promoters_shared_across_species.bed | awk '{{print $1}}')
shared_enhancers=$(wc -l {output_dir}/{tissue}_enhancers_shared_across_species.bed | awk '{{print $1}}')
total_species_1_to_species_2_promoters=$(wc -l {species_1_to_species_2_promoters} | awk '{{print $1}}')

echo "{tissue} shared promoters across species: $shared_promoters"
echo "{tissue} shared enhancers across species: $shared_enhancers"
echo "Total {species_1}-to-{species_2} conserved promoters in {tissue}: $total_species_1_to_species_2_promoters"
echo "Percentage of conserved promoters shared: $(echo "scale=2; ($shared_promoters / $total_species_1_to_species_2_promoters) * 100" | bc)%"

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    conserved_master_script: Path
    shared_master_script: Path
    output_logs: list[Path]
    error_logs: list[Path]

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
    """
    Generate scripts to run the cross-species enhancer vs promoter analysis
    
    Args:
        config: A BedtoolConfig object with TSS files and conserved regions
        
    Returns:
        A GeneratedScriptOutput object containing paths to the master script and log files
    """
    # Validate that required files exist
    assert config.species_1_tss_file is not None, "species_1_tss_file is required but not provided in config"
    assert config.species_2_tss_file is not None, "species_2_tss_file is required but not provided in config"
    assert config.species_1_to_species_2_organ_1_conserved is not None, "species_1_to_species_2_organ_1_conserved is required but not provided in config"
    assert config.species_1_to_species_2_organ_2_conserved is not None, "species_1_to_species_2_organ_2_conserved is required but not provided in config"
    assert config.species_2_to_species_1_organ_1_conserved is not None, "species_2_to_species_1_organ_1_conserved is required but not provided in config"
    assert config.species_2_to_species_1_organ_2_conserved is not None, "species_2_to_species_1_organ_2_conserved is required but not provided in config"
    
    conserved_scripts: list[Path] = []
    shared_scripts: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    
    # Define combinations for conserved regions analysis
    combinations = [
        {
            "direction": f"{config.species_1}_to_{config.species_2}",
            "tissue": config.organ_1,
            "conserved_bed": config.species_1_to_species_2_organ_1_conserved,
            "tss_file": config.species_2_tss_file  # Species 2 TSS for species_1-to-species_2
        },
        {
            "direction": f"{config.species_1}_to_{config.species_2}",
            "tissue": config.organ_2,
            "conserved_bed": config.species_1_to_species_2_organ_2_conserved,
            "tss_file": config.species_2_tss_file  # Species 2 TSS for species_1-to-species_2
        },
        {
            "direction": f"{config.species_2}_to_{config.species_1}",
            "tissue": config.organ_1,
            "conserved_bed": config.species_2_to_species_1_organ_1_conserved,
            "tss_file": config.species_1_tss_file  # Species 1 TSS for species_2-to-species_1
        },
        {
            "direction": f"{config.species_2}_to_{config.species_1}",
            "tissue": config.organ_2,
            "conserved_bed": config.species_2_to_species_1_organ_2_conserved,
            "tss_file": config.species_1_tss_file  # Species 1 TSS for species_2-to-species_1
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        direction = combo["direction"]
        tissue = combo["tissue"]
        conserved_bed = combo["conserved_bed"]
        tss_file = combo["tss_file"]
        
        script_path = config.temp_dir / f"cross_species_{direction}_{tissue}.job"
        error_log = config.output_dir / f"cross_species_{direction}_{tissue}.err.txt"
        output_log = config.output_dir / f"cross_species_{direction}_{tissue}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                species=f"{config.species_1}-{config.species_2}" if direction == f"{config.species_1}_to_{config.species_2}" else f"{config.species_2}-{config.species_1}",
                tissue=tissue,
                direction=direction,
                species_1=config.species_1,
                species_2=config.species_2,
                conserved_bed=conserved_bed,
                tss_file=tss_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        conserved_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create scripts for finding shared regions across species
    shared_combinations = [
        {
            "tissue": config.organ_1,
            "species_1_to_species_2_promoters": f"{config.output_dir}/{config.species_1}_to_{config.species_2}_{config.organ_1}_conserved_promoters.bed",
            "species_1_to_species_2_enhancers": f"{config.output_dir}/{config.species_1}_to_{config.species_2}_{config.organ_1}_conserved_enhancers.bed",
            "species_2_peak_file": config.species_2_organ_1_peak_file,
            "species_2_tss_file": config.species_2_tss_file
        },
        {
            "tissue": config.organ_2,
            "species_1_to_species_2_promoters": f"{config.output_dir}/{config.species_1}_to_{config.species_2}_{config.organ_2}_conserved_promoters.bed",
            "species_1_to_species_2_enhancers": f"{config.output_dir}/{config.species_1}_to_{config.species_2}_{config.organ_2}_conserved_enhancers.bed",
            "species_2_peak_file": config.species_2_organ_2_peak_file,
            "species_2_tss_file": config.species_2_tss_file
        }
    ]
    
    for combo in shared_combinations:
        tissue = combo["tissue"]
        species_1_to_species_2_promoters = combo["species_1_to_species_2_promoters"]
        species_1_to_species_2_enhancers = combo["species_1_to_species_2_enhancers"]
        species_2_peak_file = combo["species_2_peak_file"]
        species_2_tss_file = combo["species_2_tss_file"]
        
        script_path = config.temp_dir / f"cross_species_shared_{tissue}.job"
        error_log = config.output_dir / f"cross_species_shared_{tissue}.err.txt"
        output_log = config.output_dir / f"cross_species_shared_{tissue}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(shared_script_template.format(
                tissue=tissue,
                species_1=config.species_1,
                species_2=config.species_2,
                species_1_to_species_2_promoters=species_1_to_species_2_promoters,
                species_1_to_species_2_enhancers=species_1_to_species_2_enhancers,
                species_2_peak_file=species_2_peak_file,
                species_2_tss_file=species_2_tss_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        shared_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create master script to submit all conserved region jobs
    conserved_master_script = config.temp_dir / "submit_all_cross_species_conserved_jobs.sh"
    with open(conserved_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all cross-species conserved region classification jobs...'\n")
        for script in conserved_scripts:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    conserved_master_script.chmod(0o755)  # Make the master script executable

    # Create master script to submit all shared region jobs
    shared_master_script = config.temp_dir / "submit_all_cross_species_shared_jobs.sh"
    with open(shared_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all cross-species shared region jobs...'\n")
        for script in shared_scripts:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    shared_master_script.chmod(0o755)  # Make the master script executable

    return GeneratedScriptOutput(
        conserved_master_script=conserved_master_script,
        shared_master_script=shared_master_script,
        output_logs=output_logs,
        error_logs=error_logs,
    )

def extract_conserved_region_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract conserved region counts from output logs and save to a CSV file
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Direction", "Tissue", "Total_Regions", "Promoters", "Promoters_Pct", "Enhancers", "Enhancers_Pct"]
    
    with open(output_csv, 'w') as f:
        f.write("Direction,Tissue,Total_Regions,Promoters,Promoters_Pct,Enhancers,Enhancers_Pct\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            # Skip shared logs for this summary
            if "shared" in log_file.stem:
                continue
                
            parts = log_file.stem.replace("cross_species_", "").split("_")
            if len(parts) < 3:  # Need at least direction_species_tissue
                continue
                
            # The direction part could be species1_to_species2
            direction = "_".join(parts[0:3])  # This captures "species1_to_species2"
            tissue = parts[3]  # The tissue is after the direction
            
            total_regions = 0
            promoters = 0
            enhancers = 0
            
            with open(log_file, 'r') as log:
                for line in log:
                    if f"Total {direction} {tissue} conserved regions:" in line:
                        total_regions = int(line.strip().split()[-1])
                    elif "Conserved promoters:" in line:
                        promoters = int(line.strip().split()[-1])
                    elif "Conserved enhancers:" in line:
                        enhancers = int(line.strip().split()[-1])
            
            # Calculate percentages
            promoters_pct = round(promoters / total_regions * 100, 2) if total_regions > 0 else 0
            enhancers_pct = round(enhancers / total_regions * 100, 2) if total_regions > 0 else 0
            
            f.write(f"{direction},{tissue},{total_regions},{promoters},{promoters_pct},{enhancers},{enhancers_pct}\n")
            data.append([direction, tissue, total_regions, promoters, f"{promoters_pct}%", enhancers, f"{enhancers_pct}%"])
    
    print(f"Cross-species conserved regions counts summary saved to {output_csv}")
    print("Cross-Species Conserved Regions Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def extract_shared_counts(output_logs: list[Path], output_csv: Path, conserved_csv: Path) -> None:
    """
    Extract shared region counts from output logs and save to a CSV file
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
        conserved_csv: Path to the CSV with conserved region counts
    """
    # First load conserved region counts to calculate percentages
    conserved_counts = {}
    if conserved_csv.exists():
        with open(conserved_csv, 'r') as f:
            # Skip header
            f.readline()
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 5:  # Make sure we have enough columns
                    direction = parts[0]
                    tissue = parts[1]
                    promoters = int(parts[3])
                    enhancers = int(parts[5])
                    conserved_counts[(direction, tissue)] = (promoters, enhancers)
    
    data = []
    headers = ["Tissue", "Shared_Promoters", "Shared_Enhancers", "Promoters_Pct", "Enhancers_Pct"]
    
    with open(output_csv, 'w') as f:
        f.write("Tissue,Shared_Promoters,Shared_Enhancers,Promoters_Pct,Enhancers_Pct\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            # Only process shared region logs
            if "shared" not in log_file.name:
                continue
            
            # Parse file name - format is "cross_species_shared_TISSUE.out.txt"
            parts = log_file.stem.split('_')
            if len(parts) < 4:
                print(f"Warning: Unexpected log filename format: {log_file}, skipping")
                continue
            
            tissue = parts[3]  # Extract tissue (the last part)
            
            shared_promoters = 0
            shared_enhancers = 0
            total_species1_to_species2_promoters = 0
            species1 = ""
            species2 = ""
            
            # Parse the log file for the counts
            with open(log_file, 'r') as log:
                for line in log:
                    if "shared promoters across species:" in line:
                        shared_promoters = int(line.strip().split()[-1])
                    elif "shared enhancers across species:" in line:
                        shared_enhancers = int(line.strip().split()[-1])
                    elif "Total" in line and "conserved promoters in" in line:
                        # Line format: "Total Human-to-Mouse conserved promoters in Liver: 123"
                        parts = line.strip().split()
                        species_parts = parts[1].split("-to-")
                        if len(species_parts) == 2:
                            species1 = species_parts[0]
                            species2 = species_parts[1]
                        total_species1_to_species2_promoters = int(parts[-1])
            
            # Calculate percentages using the total conserved promoters
            promoters_pct = round((shared_promoters / total_species1_to_species2_promoters) * 100, 2) if total_species1_to_species2_promoters > 0 else 0
            
            # Get total conserved enhancers for percentage calculation
            enhancers_pct = 0
            if species1 and species2:
                direction = f"{species1}_to_{species2}"
                if (direction, tissue) in conserved_counts:
                    _, total_enhancers = conserved_counts[(direction, tissue)]
                    if total_enhancers > 0:
                        enhancers_pct = round((shared_enhancers / total_enhancers) * 100, 2)
            
            f.write(f"{tissue},{shared_promoters},{shared_enhancers},{promoters_pct},{enhancers_pct}\n")
            data.append([tissue, shared_promoters, shared_enhancers, f"{promoters_pct}%", f"{enhancers_pct}%"])
    
    print(f"Cross-species shared regions counts summary saved to {output_csv}")
    print("Cross-Species Shared Regions Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_species_enhancer_promoter_pipeline(config_path: Path) -> bool:
    """
    Run the cross-species enhancer vs promoter classification pipeline.
    
    Args:
        config_path: Path to the configuration file.
        
    Returns:
        True if the pipeline executed successfully, False otherwise.
    """
    config = load_bedtool_config(config_path, "cross_species_enhancers_vs_promoters_output_dir")
    
    script_output = generate_script(config)
    
    # Clean old output and error logs
    old_log_count = 0
    for log in script_output.output_logs + script_output.error_logs:
        if log.exists():
            log.unlink()
            old_log_count += 1
    if old_log_count > 0:
        print(f"Deleted {old_log_count} old log files")
    
    # Split log lists for monitoring
    conserved_output_logs = script_output.output_logs[:4]  # First 4 are conserved region analyses
    conserved_error_logs = script_output.error_logs[:4]
    shared_output_logs = script_output.output_logs[4:]     # Last 2 are shared region analyses
    shared_error_logs = script_output.error_logs[4:]
    
    # Define CSV output paths
    conserved_csv_output = config.output_dir / f"cross_species_conserved_regions_summary.csv"
    shared_csv_output = config.output_dir / f"cross_species_shared_regions_summary.csv"
    
    # Phase 1: Submit conserved region jobs
    print("Phase 1: Submitting conserved region jobs...")
    result = subprocess.run(["bash", str(script_output.conserved_master_script)],
                            check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False
    
    # Monitor Phase 1 jobs
    try:
        phase1_success = monitor_jobs(conserved_output_logs, conserved_error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Conserved region jobs will still run.")
        raise KeyboardInterrupt
    
    # Only proceed to Phase 2 if Phase 1 was successful
    if not phase1_success:
        print("Conserved region jobs failed. Skipping shared region analysis jobs.")
        return False
    
    # Create conserved region counts summary
    if phase1_success:
        extract_conserved_region_counts(conserved_output_logs, conserved_csv_output)
    
    # Phase 2: Submit shared region jobs
    print("\nPhase 2: Submitting shared region analysis jobs...")
    result = subprocess.run(["bash", str(script_output.shared_master_script)],
                            check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False

    # Monitor Phase 2 jobs
    try:
        print("Waiting for shared region analysis jobs to complete...")
        phase2_success = monitor_jobs(shared_output_logs, shared_error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Shared region analysis jobs will still run.")
        raise KeyboardInterrupt
    
    # If jobs completed successfully, create summary CSVs
    if phase2_success: 
        # Create shared region counts summary with percentages
        extract_shared_counts(shared_output_logs, shared_csv_output, conserved_csv_output)
    
    return phase1_success and phase2_success
