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

echo "Processing {species}: comparing {tissue1_file} and {tissue2_file}"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Find shared regions - tissue1 peaks also in tissue2
shared_file={output_dir}/{species}_{tissue1}_peaks_shared_with_{tissue2}.bed
echo "[STEP 1] Finding shared peaks for {tissue1}: $shared_file"
bedtools intersect -a {tissue1_file} -b {tissue2_file} -u > $shared_file

# Step 2: Find tissue-specific peaks
specific_file={output_dir}/{species}_{tissue1}_specific.bed
echo "[STEP 2] Finding {tissue1}-specific peaks: $specific_file"
bedtools intersect -a {tissue1_file} -b {tissue2_file} -v > $specific_file

# Step 3: Count peaks
echo "[STEP 3] Peak counts for {tissue1}:"
total_peaks=$(wc -l {tissue1_file} | awk '{{print $1}}')
shared_peaks=$(wc -l $shared_file | awk '{{print $1}}')
specific_peaks=$(wc -l $specific_file | awk '{{print $1}}')

echo "Total {tissue1} peaks: $total_peaks"
echo "Shared with {tissue2}: $shared_peaks"
echo "Specific to {tissue1}: $specific_peaks"

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    script: Path
    output_logs: list[Path]
    error_logs: list[Path]

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
    """
    Generate scripts to run the bedtools analysis for different species and tissue combinations.

    Args:
        config: A BedtoolConfig object.

    Returns:
        A GeneratedScriptOutput object containing paths to the master script and log files
    """
    script_paths: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    
    # Define the four combinations we want to analyze
    combinations = [
        {
            "species": f"{config.species_1}",
            "tissue1": config.organ_1,
            "tissue2": config.organ_2,
            "tissue1_file": config.species_1_organ_1_peak_file,
            "tissue2_file": config.species_1_organ_2_peak_file
        },
        {
            "species": f"{config.species_1}",
            "tissue1": config.organ_2,
            "tissue2": config.organ_1,
            "tissue1_file": config.species_1_organ_2_peak_file,
            "tissue2_file": config.species_1_organ_1_peak_file
        },
        {
            "species": f"{config.species_2}",
            "tissue1": config.organ_1,
            "tissue2": config.organ_2,
            "tissue1_file": config.species_2_organ_1_peak_file,
            "tissue2_file": config.species_2_organ_2_peak_file
        },
        {
            "species": f"{config.species_2}",
            "tissue1": config.organ_2,
            "tissue2": config.organ_1,
            "tissue1_file": config.species_2_organ_2_peak_file,
            "tissue2_file": config.species_2_organ_1_peak_file
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        species = combo["species"]
        tissue1 = combo["tissue1"]
        tissue2 = combo["tissue2"]
        tissue1_file = combo["tissue1_file"]
        tissue2_file = combo["tissue2_file"]
        
        script_path = config.temp_dir / f"bedtool_{species}_{tissue1}_vs_{tissue2}.job"
        error_log = config.output_dir / f"bedtool_{species}_{tissue1}_vs_{tissue2}.err.txt"
        output_log = config.output_dir / f"bedtool_{species}_{tissue1}_vs_{tissue2}.out.txt"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                species=species,
                tissue1=tissue1,
                tissue2=tissue2,
                tissue1_file=tissue1_file,
                tissue2_file=tissue2_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_tissue_comparison_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all tissue comparison jobs...'\n")
        for script in script_paths:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    master_script.chmod(0o755)  # Make the master script executable
    
    return GeneratedScriptOutput(
        script=master_script,
        output_logs=output_logs,
        error_logs=error_logs
    )

def extract_peak_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract peak counts from output logs and save them to a CSV file.
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Species", "Tissue", "Total_Peaks", "Shared_Peaks", "Specific_Peaks"]
    
    with open(output_csv, 'w') as f:
        f.write("Species,Tissue,Total_Peaks,Shared_Peaks,Specific_Peaks\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            parts = log_file.stem.replace("bedtool_", "").replace(".out", "").split("_vs_")
            if len(parts) != 2:
                continue
                
            species_tissue1 = parts[0].split("_")
            species = species_tissue1[0]
            tissue1 = species_tissue1[1]
            
            total_peaks = 0
            shared_peaks = 0
            specific_peaks = 0
            
            with open(log_file, 'r') as log:
                for line in log:
                    if f"Total {tissue1} peaks:" in line:
                        total_peaks = line.strip().split()[-1]
                    elif "Shared with" in line:
                        shared_peaks = line.strip().split()[-1]
                    elif f"Specific to {tissue1}:" in line:
                        specific_peaks = line.strip().split()[-1]
            
            f.write(f"{species},{tissue1},{total_peaks},{shared_peaks},{specific_peaks}\n")
            data.append([species, tissue1, total_peaks, shared_peaks, specific_peaks])
    
    print(f"Peak counts summary saved to {output_csv}")
    print("Peak Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_tissues_shared_vs_specific_pipeline(config_path: Path) -> bool:
    """
    Run the bedtools comparison pipeline for tissue-specific and shared regions.

    Args:
        config_path: Path to the configuration file.
    """
    config = load_bedtool_config(config_path, "cross_tissues_shared_vs_specific_output_dir")
    script_output = generate_script(config)
    script_path = script_output.script

    # Clean old output err logs
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
        raise KeyboardInterrupt
    
    # If jobs completed successfully, create a summary CSV
    if success:
        csv_output = config.output_dir / f"tissue_comparison_peak_counts_summary.csv"
        extract_peak_counts(script_output.output_logs, csv_output)
    
    return success
