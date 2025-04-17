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

echo "Processing {species} {tissue} enhancers vs promoters classification"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Sort peak files
echo "[STEP 1] Sorting peak file: {peak_file}"
sorted_file={output_dir}/{species}_{tissue}_peaks.sorted.bed
bedtools sort -i {peak_file} > $sorted_file

# Step 2: Annotate distance to TSS
echo "[STEP 2] Annotating TSS distance"
tss_annotated_file={output_dir}/{species}_{tissue}_TSS_annotated.bed
bedtools closest -a $sorted_file -b {tss_file} -d > $tss_annotated_file

# Step 3: Classify as promoters vs enhancers
echo "[STEP 3] Classifying enhancers and promoters"
promoters_file={output_dir}/{species}_{tissue}_promoters.bed
enhancers_file={output_dir}/{species}_{tissue}_enhancers.bed
awk '$NF <= 5000' $tss_annotated_file > $promoters_file
awk '$NF > 5000' $tss_annotated_file > $enhancers_file

# Step 4: Count peaks
echo "[STEP 4] Counting peaks"
total_peaks=$(wc -l $sorted_file | awk '{{print $1}}')
promoters=$(wc -l $promoters_file | awk '{{print $1}}')
enhancers=$(wc -l $enhancers_file | awk '{{print $1}}')

echo "Total {species} {tissue} peaks: $total_peaks"
echo "Promoters: $promoters"
echo "Enhancers: $enhancers"

echo "Job finished"
"""

shared_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:20:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {species} shared enhancers and promoters between tissues"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Wait for required files to be available before proceeding
echo "Checking if required files exist..."
max_wait_time=600  # Maximum wait time in seconds (10 minutes)
wait_interval=1    # Check every 1 seconds
elapsed_time=0

while [ $elapsed_time -lt $max_wait_time ]; do
    all_files_exist=true
    
    # Check tissue1 files
    if [ ! -f {tissue1_promoters} ]; then
        echo "Waiting for file: {tissue1_promoters} (not found)"
        all_files_exist=false
    fi
    if [ ! -f {tissue1_enhancers} ]; then
        echo "Waiting for file: {tissue1_enhancers} (not found)"
        all_files_exist=false
    fi
    
    # Check tissue2 files
    if [ ! -f {tissue2_promoters} ]; then
        echo "Waiting for file: {tissue2_promoters} (not found)"
        all_files_exist=false
    fi
    if [ ! -f {tissue2_enhancers} ]; then
        echo "Waiting for file: {tissue2_enhancers} (not found)"
        all_files_exist=false
    fi
    
    # If all files exist, proceed
    if [ "$all_files_exist" = true ]; then
        echo "All required files found. Proceeding with analysis."
        break
    fi
    
    # Wait before checking again
    echo "Waiting for files to be created... (${{elapsed_time}}s elapsed)"
    sleep $wait_interval
    elapsed_time=$((elapsed_time + wait_interval))
done

# If we've hit the maximum wait time, exit with an error
if [ $elapsed_time -ge $max_wait_time ]; then
    echo "Error: Maximum wait time exceeded. Required files not found within 10 minutes." >&2
    exit 1
fi

# Step 1: Find shared promoters between tissues
echo "[STEP 1] Finding shared promoters"
shared_promoters_file={output_dir}/{species}_promoters_shared_across_tissues.bed
bedtools intersect -a {tissue1_promoters} -b {tissue2_promoters} -u > $shared_promoters_file

# Step 2: Find shared enhancers between tissues
echo "[STEP 2] Finding shared enhancers"
shared_enhancers_file={output_dir}/{species}_enhancers_shared_across_tissues.bed
bedtools intersect -a {tissue1_enhancers} -b {tissue2_enhancers} -u > $shared_enhancers_file

# Step 3: Count shared elements
echo "[STEP 3] Counting shared elements"
promoters_t1=$(wc -l {tissue1_promoters} | awk '{{print $1}}')
promoters_t2=$(wc -l {tissue2_promoters} | awk '{{print $1}}')
enhancers_t1=$(wc -l {tissue1_enhancers} | awk '{{print $1}}')
enhancers_t2=$(wc -l {tissue2_enhancers} | awk '{{print $1}}')
shared_promoters=$(wc -l $shared_promoters_file | awk '{{print $1}}')
shared_enhancers=$(wc -l $shared_enhancers_file | awk '{{print $1}}')

echo "{species} {tissue1} promoters: $promoters_t1"
echo "{species} {tissue2} promoters: $promoters_t2"
echo "{species} shared promoters: $shared_promoters"
echo "{species} {tissue1} enhancers: $enhancers_t1"
echo "{species} {tissue2} enhancers: $enhancers_t2"
echo "{species} shared enhancers: $shared_enhancers"

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    script: Path
    output_logs: list[Path]
    error_logs: list[Path]

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
    """
    Generate scripts to run the enhancer vs promoter analysis
    
    Args:
        config: A BedtoolConfig object with TSS files
        
    Returns:
        A GeneratedScriptOutput object containing paths to the master script and log files
    """
    # Validate that TSS files exist
    assert config.species_1_tss_file is not None, "species_1_tss_file is required but not provided in config"
    assert config.species_2_tss_file is not None, "species_2_tss_file is required but not provided in config"
    
    script_paths: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    
    # Define the four combinations for initial classification
    combinations = [
        {
            "species": f"{config.species_1}",
            "tissue": config.organ_1, 
            "peak_file": config.species_1_organ_1_peak_file,
            "tss_file": config.species_1_tss_file
        },
        {
            "species": f"{config.species_1}",
            "tissue": config.organ_2,
            "peak_file": config.species_1_organ_2_peak_file,
            "tss_file": config.species_1_tss_file
        },
        {
            "species": f"{config.species_2}",
            "tissue": config.organ_1,
            "peak_file": config.species_2_organ_1_peak_file,
            "tss_file": config.species_2_tss_file
        },
        {
            "species": f"{config.species_2}",
            "tissue": config.organ_2,
            "peak_file": config.species_2_organ_2_peak_file,
            "tss_file": config.species_2_tss_file
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        species = combo["species"]
        tissue = combo["tissue"]
        peak_file = combo["peak_file"]
        tss_file = combo["tss_file"]
        
        script_path = config.temp_dir / f"enhancer_promoter_{species}_{tissue}.job"
        error_log = config.output_dir / f"enhancer_promoter_{species}_{tissue}.err.txt"
        output_log = config.output_dir / f"enhancer_promoter_{species}_{tissue}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                species=species,
                tissue=tissue,
                peak_file=peak_file,
                tss_file=tss_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create scripts for finding shared regions between tissues (per species)
    shared_combinations = [
        {
            "species": f"{config.species_1}",
            "tissue1": config.organ_1,
            "tissue2": config.organ_2,
        },
        {
            "species": f"{config.species_2}",
            "tissue1": config.organ_1,
            "tissue2": config.organ_2,
        }
    ]
    
    for combo in shared_combinations:
        species = combo["species"]
        tissue1 = combo["tissue1"]
        tissue2 = combo["tissue2"]
        
        script_path = config.temp_dir / f"enhancer_promoter_{species}_shared.job"
        error_log = config.output_dir / f"enhancer_promoter_{species}_shared.err.txt"
        output_log = config.output_dir / f"enhancer_promoter_{species}_shared.out.txt"
        
        with open(script_path, "w") as f:
            f.write(shared_script_template.format(
                species=species,
                tissue1=tissue1,
                tissue2=tissue2,
                tissue1_promoters=f"{config.output_dir}/{species}_{tissue1}_promoters.bed",
                tissue2_promoters=f"{config.output_dir}/{species}_{tissue2}_promoters.bed",
                tissue1_enhancers=f"{config.output_dir}/{species}_{tissue1}_enhancers.bed",
                tissue2_enhancers=f"{config.output_dir}/{species}_{tissue2}_enhancers.bed",
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_enhancer_promoter_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all enhancer-promoter classification jobs...'\n")
        for script in script_paths:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    master_script.chmod(0o755)  # Make the master script executable
    
    return GeneratedScriptOutput(
        script=master_script,
        output_logs=output_logs,
        error_logs=error_logs
    )

def extract_enhancer_promoter_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract enhancer and promoter counts from output logs and save to a CSV file
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Species", "Tissue", "Total_Peaks", "Promoters", "Enhancers"]
    
    with open(output_csv, 'w') as f:
        f.write("Species,Tissue,Total_Peaks,Promoters,Enhancers\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            # Skip shared tissue logs for this summary
            if "shared" in log_file.stem:
                continue
                
            parts = log_file.stem.replace("enhancer_promoter_", "").replace(".out", "").split("_")
            if len(parts) != 2:
                continue
                
            species = parts[0]
            tissue = parts[1]
            
            total_peaks = 0
            promoters = 0
            enhancers = 0
            
            with open(log_file, 'r') as log:
                for line in log:
                    if f"Total {species} {tissue} peaks:" in line:
                        total_peaks = line.strip().split()[-1]
                    elif "Promoters:" in line:
                        promoters = line.strip().split()[-1]
                    elif "Enhancers:" in line:
                        enhancers = line.strip().split()[-1]
            
            f.write(f"{species},{tissue},{total_peaks},{promoters},{enhancers}\n")
            data.append([species, tissue, total_peaks, promoters, enhancers])
    
    print(f"Enhancer-Promoter counts summary saved to {output_csv}")
    print("Enhancer-Promoter Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def extract_shared_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract shared enhancer and promoter counts from output logs and save to a CSV file
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Species", "Tissue1_Promoters", "Tissue2_Promoters", "Shared_Promoters", 
               "Tissue1_Enhancers", "Tissue2_Enhancers", "Shared_Enhancers"]
    
    with open(output_csv, 'w') as f:
        f.write("Species,Tissue1_Promoters,Tissue2_Promoters,Shared_Promoters,Tissue1_Enhancers,Tissue2_Enhancers,Shared_Enhancers\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            # Only process shared tissue logs
            if "shared" not in log_file.stem:
                continue
                
            parts = log_file.stem.replace("enhancer_promoter_", "").replace("_shared", "").replace(".out", "")
            species = parts
            
            tissue1_promoters = 0
            tissue2_promoters = 0
            shared_promoters = 0
            tissue1_enhancers = 0
            tissue2_enhancers = 0
            shared_enhancers = 0
            
            tissue1 = ""
            tissue2 = ""
            
            with open(log_file, 'r') as log:
                for line in log:
                    if f"{species}" in line and "promoters:" in line:
                        if "shared promoters:" in line:
                            shared_promoters = line.strip().split()[-1]
                        elif tissue1 == "":
                            parts = line.split()
                            tissue1 = parts[1]
                            tissue1_promoters = parts[-1]
                        else:
                            parts = line.split()
                            tissue2 = parts[1]
                            tissue2_promoters = parts[-1]
                    elif f"{species}" in line and "enhancers:" in line:
                        if "shared enhancers:" in line:
                            shared_enhancers = line.strip().split()[-1]
                        elif tissue1 != "" and tissue1_enhancers == 0:
                            tissue1_enhancers = line.strip().split()[-1]
                        else:
                            tissue2_enhancers = line.strip().split()[-1]
            
            f.write(f"{species},{tissue1_promoters},{tissue2_promoters},{shared_promoters},{tissue1_enhancers},{tissue2_enhancers},{shared_enhancers}\n")
            data.append([species, tissue1_promoters, tissue2_promoters, shared_promoters, 
                        tissue1_enhancers, tissue2_enhancers, shared_enhancers])
    
    print(f"Shared enhancer-promoter counts summary saved to {output_csv}")
    print("Shared Enhancer-Promoter Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_tissues_enhancer_promoter_pipeline(config_path: Path) -> bool:
    """
    Run the enhancer vs promoter classification pipeline.
    
    Args:
        config_path: Path to the configuration file.
        
    Returns:
        True if the pipeline executed successfully, False otherwise.
    """
    config = load_bedtool_config(config_path, "cross_tissues_enhancers_vs_promoters_output_dir")
    
    script_output = generate_script(config)
    script_path = script_output.script
    
    # Clean old output and error logs
    old_log_count = 0
    for log in script_output.output_logs + script_output.error_logs:
        if log.exists():
            log.unlink()
            old_log_count += 1
    if old_log_count > 0:
        print(f"Deleted {old_log_count} old log files")
    
    # Submit the jobs
    result = subprocess.run(["bash", str(script_path)], check=True, capture_output=True, text=True)
    
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False
    
    # Monitor the submitted jobs
    try:
        success = monitor_jobs(script_output.output_logs, script_output.error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Jobs will still run.")
        success = False
    
    # If jobs completed successfully, create summary CSVs
    if success:
        # Create enhancer-promoter counts summary
        csv_output = config.output_dir / f"enhancer_promoter_counts_summary.csv"
        extract_enhancer_promoter_counts(script_output.output_logs, csv_output)
        
        # Create shared counts summary
        shared_csv_output = config.output_dir / f"shared_enhancer_promoter_counts_summary.csv"
        extract_shared_counts(script_output.output_logs, shared_csv_output)
    
    return success
