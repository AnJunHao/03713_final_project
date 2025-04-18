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

echo "Processing {species_from} to {species_to} {tissue} conserved regions classification"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Annotate distance to TSS
echo "[STEP 1] Annotating TSS distance"
tss_annotated_file={output_dir}/{species_from}_to_{species_to}_{tissue}_conserved.TSS.bed
bedtools closest -a {conserved_file} \
                 -b {tss_file} \
                 -d -t first > $tss_annotated_file

# Step 2: Classify as promoters vs enhancers
echo "[STEP 2] Classifying enhancers and promoters"
promoters_file={output_dir}/{species_from}_to_{species_to}_{tissue}_conserved_promoters.bed
enhancers_file={output_dir}/{species_from}_to_{species_to}_{tissue}_conserved_enhancers.bed
awk '$NF <= 5000' $tss_annotated_file > $promoters_file
awk '$NF > 5000' $tss_annotated_file > $enhancers_file

# Step 3: Count peaks
echo "[STEP 3] Counting peaks"
total_conserved=$(wc -l {conserved_file} | awk '{{print $1}}')
promoters=$(wc -l $promoters_file | awk '{{print $1}}')
enhancers=$(wc -l $enhancers_file | awk '{{print $1}}')

echo "Total {species_from} to {species_to} {tissue} conserved regions: $total_conserved"
echo "Promoters: $promoters"
echo "Enhancers: $enhancers"

echo "Job finished"
"""

shared_regions_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {tissue} shared enhancers/promoters across species - {from_species} to {to_species}"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Find shared promoters/enhancers between species
echo "[STEP 1] Finding shared promoters/enhancers across species"
shared_promoters_file={output_dir}/{tissue}_promoters_shared_across_species_{from_species}_to_{to_species}.bed
shared_enhancers_file={output_dir}/{tissue}_enhancers_shared_across_species_{from_species}_to_{to_species}.bed

bedtools intersect -a {from_to_conserved_promoters} \
                   -b {to_native_promoters} \
                   -u -wa > $shared_promoters_file

bedtools intersect -a {from_to_conserved_enhancers} \
                   -b {to_native_enhancers} \
                   -u -wa > $shared_enhancers_file

# Step 2: Count shared elements
echo "[STEP 2] Counting shared elements"
shared_promoters=$(wc -l $shared_promoters_file | awk '{{print $1}}')
shared_enhancers=$(wc -l $shared_enhancers_file | awk '{{print $1}}')

echo "{tissue} shared promoters across species ({from_species} to {to_species}): $shared_promoters"
echo "{tissue} shared enhancers across species ({from_species} to {to_species}): $shared_enhancers"

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    classification_master_script: Path
    shared_master_script: Path
    output_logs: list[Path]
    error_logs: list[Path]

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
    """
    Generate scripts to run the cross-species enhancer vs promoter analysis
    
    Args:
        config: A BedtoolConfig object with conserved files and TSS files
        
    Returns:
        A GeneratedScriptOutput object containing paths to the master script and log files
    """
    # Validate that conserved files and TSS files exist
    assert config.species_1_tss_file is not None, "species_1_tss_file is required but not provided in config"
    assert config.species_2_tss_file is not None, "species_2_tss_file is required but not provided in config"
    assert config.species_1_to_species_2_organ_1_conserved is not None, "species_1_to_species_2_organ_1_conserved is required but not provided in config"
    assert config.species_1_to_species_2_organ_2_conserved is not None, "species_1_to_species_2_organ_2_conserved is required but not provided in config"
    assert config.species_2_to_species_1_organ_1_conserved is not None, "species_2_to_species_1_organ_1_conserved is required but not provided in config"
    assert config.species_2_to_species_1_organ_2_conserved is not None, "species_2_to_species_1_organ_2_conserved is required but not provided in config"
    
    classification_scripts: list[Path] = []
    shared_scripts: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    
    # Define the four combinations for initial classification (all directions)
    combinations = [
        {
            "species_from": config.species_1,
            "species_to": config.species_2,
            "tissue": config.organ_1,
            "conserved_file": config.species_1_to_species_2_organ_1_conserved,
            "tss_file": config.species_2_tss_file,  # Use the target species TSS file
        },
        {
            "species_from": config.species_1,
            "species_to": config.species_2,
            "tissue": config.organ_2,
            "conserved_file": config.species_1_to_species_2_organ_2_conserved,
            "tss_file": config.species_2_tss_file,
        },
        {
            "species_from": config.species_2,
            "species_to": config.species_1,
            "tissue": config.organ_1,
            "conserved_file": config.species_2_to_species_1_organ_1_conserved,
            "tss_file": config.species_1_tss_file,
        },
        {
            "species_from": config.species_2,
            "species_to": config.species_1,
            "tissue": config.organ_2,
            "conserved_file": config.species_2_to_species_1_organ_2_conserved,
            "tss_file": config.species_1_tss_file,
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        species_from = combo["species_from"]
        species_to = combo["species_to"]
        tissue = combo["tissue"]
        conserved_file = combo["conserved_file"]
        tss_file = combo["tss_file"]
        
        script_path = config.temp_dir / f"cross_species_enhancer_promoter_{species_from}_to_{species_to}_{tissue}.job"
        error_log = config.output_dir / f"cross_species_enhancer_promoter_{species_from}_to_{species_to}_{tissue}.err.txt"
        output_log = config.output_dir / f"cross_species_enhancer_promoter_{species_from}_to_{species_to}_{tissue}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                species_from=species_from,
                species_to=species_to,
                tissue=tissue,
                conserved_file=conserved_file,
                tss_file=tss_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        classification_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create scripts for finding shared regions across species (per tissue)
    shared_combinations = [
        {
            "tissue": config.organ_1,
            "from_species": config.species_1,
            "to_species": config.species_2,
            "direction": f"{config.species_1}_to_{config.species_2}"
        },
        {
            "tissue": config.organ_2,
            "from_species": config.species_1,
            "to_species": config.species_2,
            "direction": f"{config.species_1}_to_{config.species_2}"
        },
        {
            "tissue": config.organ_1,
            "from_species": config.species_2,
            "to_species": config.species_1,
            "direction": f"{config.species_2}_to_{config.species_1}"
        },
        {
            "tissue": config.organ_2,
            "from_species": config.species_2,
            "to_species": config.species_1,
            "direction": f"{config.species_2}_to_{config.species_1}"
        }
    ]
    
    for combo in shared_combinations:
        tissue = combo["tissue"]
        from_species = combo["from_species"]
        to_species = combo["to_species"]
        direction = combo["direction"]
        
        script_path = config.temp_dir / f"cross_species_shared_enhancer_promoter_{direction}_{tissue}.job"
        error_log = config.output_dir / f"cross_species_shared_enhancer_promoter_{direction}_{tissue}.err.txt"
        output_log = config.output_dir / f"cross_species_shared_enhancer_promoter_{direction}_{tissue}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(shared_regions_script_template.format(
                tissue=tissue,
                from_species=from_species,
                to_species=to_species,
                from_to_conserved_promoters=f"{config.output_dir}/{from_species}_to_{to_species}_{tissue}_conserved_promoters.bed",
                to_native_promoters=f"{config.output_dir}/{to_species}_{tissue}_promoters.bed",
                from_to_conserved_enhancers=f"{config.output_dir}/{from_species}_to_{to_species}_{tissue}_conserved_enhancers.bed",
                to_native_enhancers=f"{config.output_dir}/{to_species}_{tissue}_enhancers.bed",
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        shared_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create master script to submit all classification jobs
    classification_master_script = config.temp_dir / "submit_all_cross_species_classification_jobs.sh"
    with open(classification_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all cross-species enhancer-promoter classification jobs...'\n")
        for script in classification_scripts:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    classification_master_script.chmod(0o755)  # Make the master script executable

    # Create master script to submit all shared region jobs
    shared_master_script = config.temp_dir / "submit_all_cross_species_shared_region_jobs.sh"
    with open(shared_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all cross-species shared enhancer-promoter jobs...'\n")
        for script in shared_scripts:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    shared_master_script.chmod(0o755)  # Make the master script executable

    return GeneratedScriptOutput(
        classification_master_script=classification_master_script,
        shared_master_script=shared_master_script,
        output_logs=output_logs,
        error_logs=error_logs,
    )

def extract_classification_counts(output_logs: list[Path], output_csv: Path) -> list:
    """
    Extract enhancer and promoter counts from cross-species classification output logs and save to a CSV file
    
    Args:
        output_logs: List of output log paths
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Species From", "Species To", "Tissue", "Total Conserved", "Promoters", "Promoters %", "Enhancers", "Enhancers %"]
    
    for log_file in output_logs:
        # Only process logs from classification jobs (not shared region jobs)
        if "cross_species_enhancer_promoter" in log_file.name and "shared" not in log_file.name:
            species_from = None
            species_to = None
            tissue = None
            total_conserved = None
            promoters = None
            enhancers = None
            
            with open(log_file, "r") as f:
                for line in f:
                    if "Processing" in line and "conserved regions classification" in line:
                        parts = line.strip().split()
                        species_from = parts[1]
                        species_to = parts[3]
                        tissue = parts[4]
                    elif "Total" in line and "conserved regions" in line:
                        total_conserved = line.strip().split(":")[-1].strip()
                    elif "Promoters:" in line:
                        promoters = line.strip().split(":")[-1].strip()
                    elif "Enhancers:" in line:
                        enhancers = line.strip().split(":")[-1].strip()
            
            if species_from and species_to and tissue and total_conserved and promoters and enhancers:
                total_conserved_int = int(total_conserved)
                promoters_int = int(promoters)
                enhancers_int = int(enhancers)
                
                # Calculate percentages
                promoters_pct = round(promoters_int / total_conserved_int * 100, 2) if total_conserved_int > 0 else 0
                enhancers_pct = round(enhancers_int / total_conserved_int * 100, 2) if total_conserved_int > 0 else 0
                
                data.append([
                    species_from, 
                    species_to, 
                    tissue, 
                    total_conserved, 
                    promoters, 
                    f"{promoters_pct}%", 
                    enhancers, 
                    f"{enhancers_pct}%"
                ])
    
    # Write to CSV
    with open(output_csv, "w") as f:
        f.write(",".join(headers) + "\n")
        for row in data:
            f.write(",".join(str(item) for item in row) + "\n")
    
    # Print a summary table
    print("\nCross-Species Enhancer-Promoter Classification Results:")
    print(tabulate(data, headers=headers, tablefmt="grid"))
    
    # Return the data for use in other functions
    return data

def extract_shared_counts(output_logs: list[Path], output_csv: Path, classification_csv: Path) -> None:
    """
    Extract shared enhancer and promoter counts across species from output logs and save to a CSV file
    
    Args:
        output_logs: List of output log paths
        output_csv: Path to save the CSV file
        classification_csv: Path to the classification CSV file (for reference)
    """
    # First load the classification data to get the total counts
    classification_data = {}
    if classification_csv.exists():
        with open(classification_csv, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if i > 0:  # Skip header line
                    parts = line.strip().split(",")
                    if len(parts) >= 7:
                        species_from = parts[0]
                        species_to = parts[1]
                        tissue = parts[2]
                        promoters = int(parts[4])
                        enhancers = int(parts[6])
                        
                        key = f"{species_from}_{species_to}_{tissue}"
                        classification_data[key] = {"promoters": promoters, "enhancers": enhancers}
    
    data = []
    headers = ["Tissue", "From Species", "To Species", "Shared Promoters", "Shared Promoters %", "Shared Enhancers", "Shared Enhancers %"]
    
    for log_file in output_logs:
        # Only process logs from shared region jobs
        if "cross_species_shared_enhancer_promoter" in log_file.name:
            tissue = None
            from_species = None
            to_species = None
            shared_promoters = None
            shared_enhancers = None
            
            with open(log_file, "r") as f:
                for line in f:
                    if "Processing" in line and "shared enhancers/promoters across species" in line:
                        parts = line.strip().split()
                        tissue = parts[1]
                        # Extract from/to species from the line (format: "from_species to to_species")
                        direction_part = line.split("-")[1].strip()
                        from_species = direction_part.split("to")[0].strip()
                        to_species = direction_part.split("to")[1].strip()
                    elif "shared promoters across species" in line:
                        shared_promoters = line.strip().split(":")[-1].strip()
                    elif "shared enhancers across species" in line:
                        shared_enhancers = line.strip().split(":")[-1].strip()
            
            if tissue and from_species and to_species and shared_promoters and shared_enhancers:
                shared_promoters_int = int(shared_promoters)
                shared_enhancers_int = int(shared_enhancers)
                
                # Calculate percentages if classification data is available
                key = f"{from_species}_{to_species}_{tissue}"
                promoters_pct = "N/A"
                enhancers_pct = "N/A"
                
                if key in classification_data:
                    total_promoters = classification_data[key]["promoters"]
                    total_enhancers = classification_data[key]["enhancers"]
                    
                    if total_promoters > 0:
                        promoters_pct = f"{round(shared_promoters_int / total_promoters * 100, 2)}%"
                    if total_enhancers > 0:
                        enhancers_pct = f"{round(shared_enhancers_int / total_enhancers * 100, 2)}%"
                
                data.append([
                    tissue, 
                    from_species, 
                    to_species, 
                    shared_promoters, 
                    promoters_pct, 
                    shared_enhancers, 
                    enhancers_pct
                ])
    
    # Write to CSV
    with open(output_csv, "w") as f:
        f.write(",".join(headers) + "\n")
        for row in data:
            f.write(",".join(str(item) for item in row) + "\n")
    
    # Print a summary table
    print("\nCross-Species Shared Enhancer-Promoter Results:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_species_enhancer_promoter_pipeline(config_path: Path) -> bool:
    """
    Run the cross-species enhancer promoter pipeline
    
    Args:
        config_path: Path to the config file
        
    Returns:
        True if the pipeline completed successfully, False otherwise
    """
    # Load config file
    config = load_bedtool_config(config_path, "cross_species_enhancers_vs_promoters_output_dir")
    
    # Generate scripts
    script_output = generate_script(config)
    
    # Run classification jobs first
    subprocess.run([script_output.classification_master_script])
    
    # Wait for classification jobs to complete
    print("\nWaiting for cross-species classification jobs to complete...")
    if not monitor_jobs(script_output.output_logs[:4], script_output.error_logs[:4]):
        print("Error: Some cross-species classification jobs failed.")
        return False
    
    # Extract and save classification counts
    classification_csv = config.output_dir / "cross_species_classification_counts.csv"
    extract_classification_counts(script_output.output_logs, classification_csv)
    
    # Run shared region jobs
    subprocess.run([script_output.shared_master_script])
    
    # Wait for shared region jobs to complete
    print("\nWaiting for cross-species shared region jobs to complete...")
    if not monitor_jobs(script_output.output_logs[4:], script_output.error_logs[4:]):
        print("Error: Some cross-species shared region jobs failed.")
        return False
    
    # Extract and save shared counts
    shared_csv = config.output_dir / "cross_species_shared_counts.csv"
    extract_shared_counts(script_output.output_logs, shared_csv, classification_csv)
    
    print("\nCross-species enhancer promoter pipeline completed successfully!")
    return True
