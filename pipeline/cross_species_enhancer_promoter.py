from pathlib import Path
import subprocess
from typing import NamedTuple, Literal
from pipeline.monitor import monitor_jobs
from pipeline.bedtool_preprocess import BedtoolConfig, load_bedtool_config
from tabulate import tabulate
import yaml

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:20:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {species_from}_{organ} to {species_to}_{organ} enhancer-promoter analysis"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Find nearest TSS for conserved regions
echo "[STEP 1] Finding nearest TSS for conserved regions"
bedtools closest -a {conserved_file} \\
                -b {tss_file} \\
                -d > {output_dir}/{prefix}_conserved.TSS.bed

# Step 2: Categorize as promoters (<=5000bp from TSS) or enhancers (>5000bp from TSS)
echo "[STEP 2] Categorizing as promoters or enhancers"
awk '$NF <= 5000' {output_dir}/{prefix}_conserved.TSS.bed > {output_dir}/{prefix}_conserved_promoters.bed
awk '$NF > 5000' {output_dir}/{prefix}_conserved.TSS.bed > {output_dir}/{prefix}_conserved_enhancers.bed

# Step 3: Report counts
echo "[STEP 3] Counting results"
total_regions=$(wc -l {conserved_file} | awk '{{print $1}}')
promoters=$(wc -l {output_dir}/{prefix}_conserved_promoters.bed | awk '{{print $1}}')
enhancers=$(wc -l {output_dir}/{prefix}_conserved_enhancers.bed | awk '{{print $1}}')

echo "Total conserved regions: $total_regions"
echo "Promoters (<=5000bp from TSS): $promoters"
echo "Enhancers (>5000bp from TSS): $enhancers"

echo "Job finished"
"""

shared_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:15:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing shared {organ} enhancers/promoters across species analysis"

# Create output directory if it doesn't exist
mkdir -p {output_dir}

# Step 1: Find shared promoters across species
echo "[STEP 1] Finding shared promoters across species"
bedtools intersect -a {species_1_to_species_2_promoters} \\
                  -b {species_2_native_promoters} \\
                  -u > {output_dir}/{organ}_promoters_shared_across_species.bed

# Step 2: Find shared enhancers across species
echo "[STEP 2] Finding shared enhancers across species"
bedtools intersect -a {species_1_to_species_2_enhancers} \\
                  -b {species_2_native_enhancers} \\
                  -u > {output_dir}/{organ}_enhancers_shared_across_species.bed

# Step 3: Report counts
echo "[STEP 3] Counting results"
total_promoters=$(wc -l {species_1_to_species_2_promoters} | awk '{{print $1}}')
shared_promoters=$(wc -l {output_dir}/{organ}_promoters_shared_across_species.bed | awk '{{print $1}}')
total_enhancers=$(wc -l {species_1_to_species_2_enhancers} | awk '{{print $1}}')
shared_enhancers=$(wc -l {output_dir}/{organ}_enhancers_shared_across_species.bed | awk '{{print $1}}')

echo "Total conserved promoters: $total_promoters"
echo "Shared promoters across species: $shared_promoters"
echo "Shared promoters percentage: $(echo "scale=2; $shared_promoters/$total_promoters*100" | bc)%"
echo "Total conserved enhancers: $total_enhancers"
echo "Shared enhancers across species: $shared_enhancers"
echo "Shared enhancers percentage: $(echo "scale=2; $shared_enhancers/$total_enhancers*100" | bc)%"

echo "Job finished"
"""

class EnhancerPromoterOutput(NamedTuple):
    script: Path
    output_logs: list[Path]
    error_logs: list[Path]
    promoter_files: dict[str, Path]
    enhancer_files: dict[str, Path]

def generate_scripts(config: BedtoolConfig) -> EnhancerPromoterOutput:
    """
    Generate scripts to run the enhancer-promoter analysis for different species and organ combinations.

    Args:
        config: A BedtoolConfig object.

    Returns:
        An EnhancerPromoterOutput object containing paths to scripts and output files
    """
    script_paths: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    promoter_files: dict[str, Path] = {}
    enhancer_files: dict[str, Path] = {}
    
    # Define the four combinations for conserved regions analysis
    combinations = [
        {
            "prefix": f"{config.species_1}_to_{config.species_2}_{config.organ_1}",
            "species_from": config.species_1,
            "species_to": config.species_2,
            "organ": config.organ_1,
            "conserved_file": config.species_1_to_species_2_organ_1_conserved,
            "tss_file": config.species_2_tss_file
        },
        {
            "prefix": f"{config.species_1}_to_{config.species_2}_{config.organ_2}",
            "species_from": config.species_1,
            "species_to": config.species_2,
            "organ": config.organ_2,
            "conserved_file": config.species_1_to_species_2_organ_2_conserved,
            "tss_file": config.species_2_tss_file
        },
        {
            "prefix": f"{config.species_2}_to_{config.species_1}_{config.organ_1}",
            "species_from": config.species_2,
            "species_to": config.species_1,
            "organ": config.organ_1,
            "conserved_file": config.species_2_to_species_1_organ_1_conserved,
            "tss_file": config.species_1_tss_file
        },
        {
            "prefix": f"{config.species_2}_to_{config.species_1}_{config.organ_2}",
            "species_from": config.species_2,
            "species_to": config.species_1,
            "organ": config.organ_2,
            "conserved_file": config.species_2_to_species_1_organ_2_conserved,
            "tss_file": config.species_1_tss_file
        }
    ]
    
    # Create scripts for each combination
    for combo in combinations:
        prefix = combo["prefix"]
        
        script_path = config.temp_dir / f"enhancer_promoter_{prefix}.job"
        error_log = config.output_dir / f"enhancer_promoter_{prefix}.err.txt"
        output_log = config.output_dir / f"enhancer_promoter_{prefix}.out.txt"
        promoter_file = config.output_dir / f"{prefix}_conserved_promoters.bed"
        enhancer_file = config.output_dir / f"{prefix}_conserved_enhancers.bed"
        
        with open(script_path, "w") as f:
            f.write(script_template.format(
                prefix=prefix,
                species_from=combo["species_from"],
                species_to=combo["species_to"],
                organ=combo["organ"],
                conserved_file=combo["conserved_file"],
                tss_file=combo["tss_file"],
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
        
        # Store the output files for later use
        promoter_files[prefix] = promoter_file
        enhancer_files[prefix] = enhancer_file
    
    # Create shared analysis scripts for each organ
    for organ in [config.organ_1, config.organ_2]:
        prefix = f"shared_{organ}"
        script_path = config.temp_dir / f"enhancer_promoter_{prefix}.job"
        error_log = config.output_dir / f"enhancer_promoter_{prefix}.err.txt"
        output_log = config.output_dir / f"enhancer_promoter_{prefix}.out.txt"
        
        species_1_to_species_2_prefix = f"{config.species_1}_to_{config.species_2}_{organ}"
        species_2_native_prefix = f"{config.species_2}_{organ}"
        
        with open(script_path, "w") as f:
            f.write(shared_script_template.format(
                organ=organ,
                species_1_to_species_2_promoters=promoter_files[species_1_to_species_2_prefix],
                species_1_to_species_2_enhancers=enhancer_files[species_1_to_species_2_prefix],
                species_2_native_promoters=config.output_dir / f"{species_2_native_prefix}_promoters.bed",
                species_2_native_enhancers=config.output_dir / f"{species_2_native_prefix}_enhancers.bed",
                output_dir=config.output_dir,
                error_log=error_log,
                output_log=output_log
            ))
        
        script_paths.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Create native species enhancer/promoter scripts
    for species, species_name in [(1, config.species_1), (2, config.species_2)]:
        for organ in [config.organ_1, config.organ_2]:
            prefix = f"{species_name}_{organ}"
            script_path = config.temp_dir / f"enhancer_promoter_{prefix}.job"
            error_log = config.output_dir / f"enhancer_promoter_{prefix}.err.txt"
            output_log = config.output_dir / f"enhancer_promoter_{prefix}.out.txt"
            
            if species == 1:
                peak_file = config.species_1_organ_1_peak_file if organ == config.organ_1 else config.species_1_organ_2_peak_file
                tss_file = config.species_1_tss_file
            else:
                peak_file = config.species_2_organ_1_peak_file if organ == config.organ_1 else config.species_2_organ_2_peak_file
                tss_file = config.species_2_tss_file
            
            with open(script_path, "w") as f:
                f.write(f"""#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Processing {species_name}_{organ} native enhancer-promoter analysis"

# Create output directory if it doesn't exist
mkdir -p {config.output_dir}

# Step 1: Find nearest TSS for peaks
echo "[STEP 1] Finding nearest TSS for peaks"
bedtools closest -a {peak_file} \\
                -b {tss_file} \\
                -d > {config.output_dir}/{prefix}_TSS_annotated.bed

# Step 2: Categorize as promoters (<=5000bp from TSS) or enhancers (>5000bp from TSS)
echo "[STEP 2] Categorizing as promoters or enhancers"
awk '$NF <= 5000' {config.output_dir}/{prefix}_TSS_annotated.bed > {config.output_dir}/{prefix}_promoters.bed
awk '$NF > 5000' {config.output_dir}/{prefix}_TSS_annotated.bed > {config.output_dir}/{prefix}_enhancers.bed

# Step 3: Report counts
echo "[STEP 3] Counting results"
total_peaks=$(wc -l {peak_file} | awk '{{print $1}}')
promoters=$(wc -l {config.output_dir}/{prefix}_promoters.bed | awk '{{print $1}}')
enhancers=$(wc -l {config.output_dir}/{prefix}_enhancers.bed | awk '{{print $1}}')

echo "Total peaks: $total_peaks"
echo "Promoters (<=5000bp from TSS): $promoters"
echo "Enhancers (>5000bp from TSS): $enhancers"
echo "Promoters percentage: $(echo "scale=2; $promoters/$total_peaks*100" | bc)%"
echo "Enhancers percentage: $(echo "scale=2; $enhancers/$total_peaks*100" | bc)%"

echo "Job finished"
""")
            
            script_paths.append(script_path)
            output_logs.append(output_log)
            error_logs.append(error_log)
    
    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_enhancer_promoter_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all enhancer-promoter analysis jobs...'\n")
        # First run the native species jobs
        for script in script_paths[:8]:  # Assuming the first 8 scripts are for native species
            f.write(f"sbatch {script}\n")
        f.write("echo 'Native species jobs submitted. Waiting for completion...'\n")
        f.write("sleep 60\n")  # Wait for the native species jobs to complete
        # Then run the conserved region jobs
        for script in script_paths[8:12]:  # Assuming the next 4 scripts are for conserved regions
            f.write(f"sbatch {script}\n")
        f.write("echo 'Conserved region jobs submitted. Waiting for completion...'\n")
        f.write("sleep 60\n")  # Wait for the conserved region jobs to complete
        # Finally run the shared analysis jobs
        for script in script_paths[12:]:  # Assuming the last 2 scripts are for shared analysis
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    
    master_script.chmod(0o755)  # Make the master script executable
    
    return EnhancerPromoterOutput(
        script=master_script,
        output_logs=output_logs,
        error_logs=error_logs,
        promoter_files=promoter_files,
        enhancer_files=enhancer_files
    )

def extract_enhancer_promoter_counts(output_logs: list[Path], output_csv: Path) -> None:
    """
    Extract enhancer and promoter counts from output logs and save them to a CSV file.
    
    Args:
        output_logs: List of paths to output log files
        output_csv: Path to save the CSV file
    """
    data = []
    headers = ["Analysis", "Total_Peaks", "Promoters", "Promoters_Pct", "Enhancers", "Enhancers_Pct"]
    
    with open(output_csv, 'w') as f:
        f.write("Analysis,Total_Peaks,Promoters,Promoters_Pct,Enhancers,Enhancers_Pct\n")
        
        for log_file in output_logs:
            if not log_file.exists():
                print(f"Warning: Log file {log_file} does not exist, skipping")
                continue
                
            prefix = log_file.stem.replace("enhancer_promoter_", "").replace(".out", "")
            total_peaks = 0
            promoters = 0
            enhancers = 0
            promoters_pct = 0
            enhancers_pct = 0
            
            with open(log_file, 'r') as log:
                for line in log:
                    if "Total peaks:" in line or "Total conserved regions:" in line:
                        total_peaks = int(line.strip().split()[-1])
                    elif "Promoters (<=5000bp from TSS):" in line:
                        promoters = int(line.strip().split()[-1])
                    elif "Enhancers (>5000bp from TSS):" in line:
                        enhancers = int(line.strip().split()[-1])
                    elif "Promoters percentage:" in line:
                        promoters_pct = float(line.strip().split()[-1].strip('%'))
                    elif "Enhancers percentage:" in line:
                        enhancers_pct = float(line.strip().split()[-1].strip('%'))
                    elif "Shared promoters percentage:" in line:
                        promoters_pct = float(line.strip().split()[-1].strip('%'))
                    elif "Shared enhancers percentage:" in line:
                        enhancers_pct = float(line.strip().split()[-1].strip('%'))
            
            if promoters_pct == 0 and total_peaks > 0:
                promoters_pct = round(promoters / total_peaks * 100, 2)
            if enhancers_pct == 0 and total_peaks > 0:
                enhancers_pct = round(enhancers / total_peaks * 100, 2)
            
            f.write(f"{prefix},{total_peaks},{promoters},{promoters_pct},{enhancers},{enhancers_pct}\n")
            data.append([prefix, total_peaks, promoters, f"{promoters_pct}%", enhancers, f"{enhancers_pct}%"])
    
    print(f"Enhancer-Promoter summary saved to {output_csv}")
    print("Enhancer-Promoter Summary:")
    print(tabulate(data, headers=headers, tablefmt="grid"))

def run_cross_species_enhancer_promoter_pipeline(config_path: Path, do_not_submit: bool = False) -> bool:
    """
    Run the cross-species enhancer-promoter analysis pipeline.

    Args:
        config_path: Path to the configuration file.
        do_not_submit: If True, generate scripts but don't submit jobs.
        
    Returns:
        True if the pipeline ran successfully, False otherwise.
    """
    config = load_bedtool_config(config_path, "cross_species_enhancers_vs_promoters_output_dir")
    script_output = generate_scripts(config)
    script_path = script_output.script

    if do_not_submit:
        return True

    # Clean old output logs
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
        csv_output = config.output_dir / f"{config.species_1}_{config.species_2}_enhancer_promoter_summary.csv"
        extract_enhancer_promoter_counts(script_output.output_logs, csv_output)
    
    return success
