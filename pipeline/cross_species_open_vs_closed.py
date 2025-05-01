from pathlib import Path
import subprocess
from typing import NamedTuple, Literal
from pipeline.monitor import monitor_jobs
from pipeline.bedtool_preprocess import Config, load_bedtool_config
from tabulate import tabulate
from pipeline.utils import update_config

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {prefix}: {halper_file} vs {native_file}"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Intersect for conserved (ortholog open)
conserved_file={output_dir}/{prefix}_conserved.bed
echo "[STEP 1] Finding conserved peaks: $conserved_file"
bedtools intersect -a {halper_file} -b {native_file} -u -wa > $conserved_file

# Step 2: Intersect for non-conserved (ortholog closed)
closed_file={output_dir}/{prefix}_closed.bed
echo "[STEP 2] Finding non-conserved peaks: $closed_file"
bedtools intersect -a {halper_file} -b {native_file} -v > $closed_file

# Step 3: Count peaks
echo "[STEP 3] Peak counts:"
lifted_total=$(wc -l {halper_file} | awk '{{print $1}}')
open_total=$(wc -l $conserved_file | awk '{{print $1}}')
closed_total=$(wc -l $closed_file | awk '{{print $1}}')

echo "Total lifted peaks: $lifted_total"
echo "Open peaks: $open_total"
echo "Closed peaks: $closed_total"

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    script: Path
    output_logs: list[Path]
    error_logs: list[Path]
    config_entries: dict[str, Path]

def generate_script(config: Config) -> GeneratedScriptOutput:
    """
    Generate scripts to run the bedtools analysis for different species and organ combinations.

    Args:
        config: A BedtoolConfig object.

    Returns:
        A GeneratedScriptOutput object containing paths to the master script and log files
    """
    script_paths: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    config_entries: dict[str, Path] = {}
    
    # Define the four combinations we want to analyze
    combinations = [
        {
            "prefix": f"{config.species_1}_{config.organ_1}_to_{config.species_2}_{config.organ_1}",
            "halper_file": config.species_1_organ_1_to_species_2,
            "native_file": config.species_2_organ_1_peak_file,
            "config_entry": "species_1_to_species_2_organ_1_conserved"
        },
        {
            "prefix": f"{config.species_1}_{config.organ_2}_to_{config.species_2}_{config.organ_2}",
            "halper_file": config.species_1_organ_2_to_species_2,
            "native_file": config.species_2_organ_2_peak_file,
            "config_entry": "species_1_to_species_2_organ_2_conserved"
        },
        {
            "prefix": f"{config.species_2}_{config.organ_1}_to_{config.species_1}_{config.organ_1}",
            "halper_file": config.species_2_organ_1_to_species_1,
            "native_file": config.species_1_organ_1_peak_file,
            "config_entry": "species_2_to_species_1_organ_1_conserved"
        },
        {
            "prefix": f"{config.species_2}_{config.organ_2}_to_{config.species_1}_{config.organ_2}",
            "halper_file": config.species_2_organ_2_to_species_1,
            "native_file": config.species_1_organ_2_peak_file,
            "config_entry": "species_2_to_species_1_organ_2_conserved"
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        prefix = combo["prefix"]
        halper_file = combo["halper_file"]
        native_file = combo["native_file"]
        
        script_path = config.temp_dir / f"bedtool_{prefix}.job"
        error_log = config.output_dir / f"bedtool_{prefix}.err.txt"
        output_log = config.output_dir / f"bedtool_{prefix}.out.txt"
        conserved_file = config.output_dir / f"{prefix}_conserved.bed"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                prefix=prefix,
                halper_file=halper_file,
                native_file=native_file,
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
        config_entries[combo["config_entry"]] = conserved_file

    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_bedtool_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all bedtool jobs...'\n")
        for script in script_paths:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    master_script.chmod(0o755)  # Make the master script executable
    
    return GeneratedScriptOutput(
        script=master_script,
        output_logs=output_logs,
        error_logs=error_logs,
        config_entries=config_entries
    )

def extract_peak_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract peak counts from output logs and save them to a CSV file.
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Job", "Total_Lifted_Peaks", "Open_Peaks", "Open_Pct", "Closed_Peaks", "Closed_Pct"]
    
    with open(output_csv, 'w') as f:
        f.write("Job,Total_Lifted_Peaks,Open_Peaks,Open_Pct,Closed_Peaks,Closed_Pct\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            prefix = log_file.stem.replace("bedtool_", "").replace(".out", "")
            total_lifted = 0
            open_peaks = 0
            closed_peaks = 0
            
            with open(log_file, 'r') as log:
                for line in log:
                    if "Total lifted peaks:" in line:
                        total_lifted = int(line.strip().split()[-1])
                    elif "Open peaks:" in line:
                        open_peaks = int(line.strip().split()[-1])
                    elif "Closed peaks:" in line:
                        closed_peaks = int(line.strip().split()[-1])
            
            # Calculate percentages
            open_pct = round(open_peaks / total_lifted * 100, 2) if total_lifted > 0 else 0
            closed_pct = round(closed_peaks / total_lifted * 100, 2) if total_lifted > 0 else 0
            
            f.write(f"{prefix},{total_lifted},{open_peaks},{open_pct},{closed_peaks},{closed_pct}\n")
            data.append([prefix, total_lifted, open_peaks, f"{open_pct}%", closed_peaks, f"{closed_pct}%"])
    
    print(f"Peak counts summary saved to {output_csv}")
    print("Peak Counts Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_species_open_vs_closed_pipeline(config_path: Path, do_not_submit: bool = False) -> bool:
    """
    Run the bedtools comparison pipeline.

    Args:
        config_path: Path to the configuration file.
    """
    config = load_bedtool_config(config_path, "cross_species_open_vs_closed_output_dir")
    script_output = generate_script(config)
    script_path = script_output.script
    update_config(config_path, script_output.config_entries, backup_suffix=".backup03")

    if do_not_submit:
        return True

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
        csv_output = config.output_dir / f"{config.species_1}_{config.species_2}_peak_counts_summary.csv"
        extract_peak_counts(script_output.output_logs, csv_output)
    
    return success