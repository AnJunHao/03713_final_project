from pathlib import Path
import subprocess
from typing import NamedTuple
from pipeline.monitor import monitor_jobs
from pipeline.bedtool_preprocess import Config, load_bedtool_config
from tabulate import tabulate

faidx_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load samtools

echo "Creating genome index for {genome_fasta}"
samtools faidx {genome_fasta}

echo "Job finished"
"""

meme_chip_script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 3:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load MEME-suite/5.4.1
module load bedtools
module load samtools

echo "Processing MEME-chip analysis for {input_name}"

# Create output directory if it doesn't exist
mkdir -p {output_dir}
touch {output_dir}/{output_name}.fa

# Map the chromosome coordinates in the bed file to the corresponding nucleotide regions in the genome
echo "Running bedtools getfasta for {input_bed}"
bedtools getfasta -fi {genome_fasta} -bed {input_bed} -s -fo {output_dir}/{output_name}.fa

# Run MEME-chip on the generated fasta file
echo "Running MEME-chip on {output_name}.fa"
cd {output_dir}
meme-chip {output_name}.fa

echo "Job finished"
"""

class GeneratedMemeChipScriptOutput(NamedTuple):
    faidx_master_script: Path
    meme_chip_master_script: Path
    output_logs: list[Path]
    error_logs: list[Path]
    meme_chip_output_dirs: dict[str, Path]

def generate_meme_chip_scripts(config: Config) -> GeneratedMemeChipScriptOutput:
    """
    Generate scripts to run the MEME-chip analysis
    
    Args:
        config: A BedtoolConfig object
        
    Returns:
        A GeneratedMemeChipScriptOutput object containing paths to the master scripts and log files
    """
    # Check for genome fasta files
    assert config.species_1_genome_fasta is not None, "species_1_genome_fasta is required but not provided in config"
    assert config.species_2_genome_fasta is not None, "species_2_genome_fasta is required but not provided in config"
    
    # Check for conserved enhancer files
    required_files = [
        ("species_1_to_species_2_organ_1_conserved_enhancers", config.species_1_to_species_2_organ_1_conserved_enhancers),
        ("species_1_to_species_2_organ_2_conserved_enhancers", config.species_1_to_species_2_organ_2_conserved_enhancers),
        ("species_2_to_species_1_organ_1_conserved_enhancers", config.species_2_to_species_1_organ_1_conserved_enhancers),
        ("species_2_to_species_1_organ_2_conserved_enhancers", config.species_2_to_species_1_organ_2_conserved_enhancers),
        ("species_1_to_species_2_organ_1_conserved_promoters", config.species_1_to_species_2_organ_1_conserved_promoters),
        ("species_1_to_species_2_organ_2_conserved_promoters", config.species_1_to_species_2_organ_2_conserved_promoters),
        ("species_2_to_species_1_organ_1_conserved_promoters", config.species_2_to_species_1_organ_1_conserved_promoters),
        ("species_2_to_species_1_organ_2_conserved_promoters", config.species_2_to_species_1_organ_2_conserved_promoters),
    ]
    
    for name, file_path in required_files:
        assert file_path is not None, f"{name} is required but not provided in config"
    
    faidx_scripts: list[Path] = []
    meme_chip_scripts: list[Path] = []
    output_logs: list[Path] = []
    error_logs: list[Path] = []
    meme_chip_output_dirs: dict[str, Path] = {}
    
    # Generate faidx scripts for each genome
    faidx_jobs = [
        {
            "genome_fasta": config.species_1_genome_fasta,
            "species": config.species_1
        },
        {
            "genome_fasta": config.species_2_genome_fasta,
            "species": config.species_2
        }
    ]
    
    for job in faidx_jobs:
        genome_fasta = job["genome_fasta"]
        species = job["species"]
        script_path = config.temp_dir / f"faidx_{species}_genome.job"
        error_log = config.output_dir / f"faidx_{species}_genome.err.txt"
        output_log = config.output_dir / f"faidx_{species}_genome.out.txt"
        
        with open(script_path, "w") as f:
            f.write(faidx_script_template.format(
                genome_fasta=genome_fasta,
                error_log=error_log,
                output_log=output_log
            ))
        
        faidx_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
    
    # Generate MEME-chip scripts for conserved regions
    meme_chip_jobs = [
        {
            "input_name": f"{config.species_1}_to_{config.species_2}_{config.organ_1}_conserved_enhancers",
            "input_bed": config.species_1_to_species_2_organ_1_conserved_enhancers,
            "output_name": f"{config.species_1}_to_{config.species_2}_{config.organ_1}_enhancers",
            "genome_fasta": config.species_2_genome_fasta
        },
        {
            "input_name": f"{config.species_1}_to_{config.species_2}_{config.organ_2}_conserved_enhancers",
            "input_bed": config.species_1_to_species_2_organ_2_conserved_enhancers,
            "output_name": f"{config.species_1}_to_{config.species_2}_{config.organ_2}_enhancers",
            "genome_fasta": config.species_2_genome_fasta
        },
        {
            "input_name": f"{config.species_2}_to_{config.species_1}_{config.organ_1}_conserved_enhancers",
            "input_bed": config.species_2_to_species_1_organ_1_conserved_enhancers,
            "output_name": f"{config.species_2}_to_{config.species_1}_{config.organ_1}_enhancers",
            "genome_fasta": config.species_1_genome_fasta
        },
        {
            "input_name": f"{config.species_2}_to_{config.species_1}_{config.organ_2}_conserved_enhancers",
            "input_bed": config.species_2_to_species_1_organ_2_conserved_enhancers,
            "output_name": f"{config.species_2}_to_{config.species_1}_{config.organ_2}_enhancers",
            "genome_fasta": config.species_1_genome_fasta
        },
        {
            "input_name": f"{config.species_1}_to_{config.species_2}_{config.organ_1}_conserved_promoters",
            "input_bed": config.species_1_to_species_2_organ_1_conserved_promoters,
            "output_name": f"{config.species_1}_to_{config.species_2}_{config.organ_1}_promoters",
            "genome_fasta": config.species_2_genome_fasta
        },
        {
            "input_name": f"{config.species_1}_to_{config.species_2}_{config.organ_2}_conserved_promoters",
            "input_bed": config.species_1_to_species_2_organ_2_conserved_promoters,
            "output_name": f"{config.species_1}_to_{config.species_2}_{config.organ_2}_promoters",
            "genome_fasta": config.species_2_genome_fasta
        },
        {
            "input_name": f"{config.species_2}_to_{config.species_1}_{config.organ_1}_conserved_promoters",
            "input_bed": config.species_2_to_species_1_organ_1_conserved_promoters,
            "output_name": f"{config.species_2}_to_{config.species_1}_{config.organ_1}_promoters",
            "genome_fasta": config.species_1_genome_fasta
        },
        {
            "input_name": f"{config.species_2}_to_{config.species_1}_{config.organ_2}_conserved_promoters",
            "input_bed": config.species_2_to_species_1_organ_2_conserved_promoters,
            "output_name": f"{config.species_2}_to_{config.species_1}_{config.organ_2}_promoters",
            "genome_fasta": config.species_1_genome_fasta
        }
    ]
    
    for job in meme_chip_jobs:
        input_name = job["input_name"]
        input_bed = job["input_bed"]
        output_name = job["output_name"]
        genome_fasta = job["genome_fasta"]
        
        script_path = config.temp_dir / f"meme_chip_{output_name}.job"
        error_log = config.output_dir / f"meme_chip_{output_name}.err.txt"
        output_log = config.output_dir / f"meme_chip_{output_name}.out.txt"
        output_dir = config.output_dir / f"meme_{output_name}"
        
        with open(script_path, "w") as f:
            f.write(meme_chip_script_template.format(
                input_name=input_name,
                input_bed=input_bed,
                output_name=output_name,
                output_dir=output_dir,
                genome_fasta=genome_fasta,
                error_log=error_log,
                output_log=output_log
            ))
        
        meme_chip_scripts.append(script_path)
        output_logs.append(output_log)
        error_logs.append(error_log)
        meme_chip_output_dirs[output_name] = output_dir
    
    # Create master script for faidx jobs
    faidx_master_script = config.temp_dir / "submit_all_faidx_jobs.sh"
    with open(faidx_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting genome indexing jobs...'\n")
        for script in faidx_scripts:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All faidx jobs submitted'\n")
    
    faidx_master_script.chmod(0o755)  # Make the master script executable
    
    # Create master script for MEME-chip jobs
    meme_chip_master_script = config.temp_dir / "submit_all_meme_chip_jobs.sh"
    with open(meme_chip_master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting MEME-chip analysis jobs...'\n")
        
        # Organize scripts by enhancers/promoters for each species-organ combination
        # This helps ensure we have at most 4 parallel jobs (one worker per enhancer/promoter pair)
        job_groups = [
            meme_chip_scripts[0:2],  # First species enhancers (2 organs)
            meme_chip_scripts[2:4],  # Second species enhancers (2 organs)
            meme_chip_scripts[4:6],  # First species promoters (2 organs)
            meme_chip_scripts[6:8],  # Second species promoters (2 organs)
        ]
        
        for i, group in enumerate(job_groups):
            f.write(f"# Worker {i+1}\n")
            for script in group:
                f.write(f"sbatch {script}\n")
            f.write("\n")
        
        f.write("echo 'All MEME-chip jobs submitted'\n")
    
    meme_chip_master_script.chmod(0o755)  # Make the master script executable
    
    return GeneratedMemeChipScriptOutput(
        faidx_master_script=faidx_master_script,
        meme_chip_master_script=meme_chip_master_script,
        output_logs=output_logs,
        error_logs=error_logs,
        meme_chip_output_dirs=meme_chip_output_dirs
    )

def run_meme_chip_analysis_pipeline(config_path: Path) -> bool:
    """
    Run the MEME-chip analysis pipeline.
    
    Args:
        config_path: Path to the configuration file.
        
    Returns:
        True if the pipeline executed successfully, False otherwise.
    """
    config = load_bedtool_config(config_path, "meme_chip_output_dir")
    
    script_output = generate_meme_chip_scripts(config)
    
    # Clean old output and error logs
    old_log_count = 0
    for log in script_output.output_logs + script_output.error_logs:
        if log.exists():
            log.unlink()
            old_log_count += 1
    if old_log_count > 0:
        print(f"Deleted {old_log_count} old log files")
    
    # Split log lists for monitoring
    faidx_output_logs = script_output.output_logs[:2]  # First 2 are faidx jobs
    faidx_error_logs = script_output.error_logs[:2]
    meme_chip_output_logs = script_output.output_logs[2:]  # Rest are MEME-chip jobs
    meme_chip_error_logs = script_output.error_logs[2:]
    
    # Phase 1: Submit faidx jobs
    print("Phase 1: Submitting genome indexing jobs...")
    result = subprocess.run(["bash", str(script_output.faidx_master_script)],
                            check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False
    
    # Monitor Phase 1 jobs
    try:
        phase1_success = monitor_jobs(faidx_output_logs, faidx_error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Genome indexing jobs will still run.")
        raise KeyboardInterrupt
    
    # Only proceed to Phase 2 if Phase 1 was successful
    if not phase1_success:
        print("Genome indexing jobs failed. Skipping MEME-chip analysis jobs.")
        return False
    
    # Phase 2: Submit MEME-chip jobs
    print("\nPhase 2: Submitting MEME-chip analysis jobs...")
    result = subprocess.run(["bash", str(script_output.meme_chip_master_script)],
                            check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False

    # Monitor Phase 2 jobs
    try:
        print("Waiting for MEME-chip analysis jobs to complete...")
        phase2_success = monitor_jobs(meme_chip_output_logs, meme_chip_error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. MEME-chip analysis jobs will still run.")
        raise KeyboardInterrupt
    
    # Print summary of MEME-chip analysis
    if phase2_success:
        print("\nMEME-chip Analysis Summary:")
        print("=" * 50)
        print("Completed MEME-chip analysis for the following files:")
        
        data = []
        headers = ["Species Direction", "Organ", "Type", "Output Directory"]
        
        for output_name, output_dir in script_output.meme_chip_output_dirs.items():
            parts = output_name.split("_")
            species_direction = f"{parts[0]} to {parts[2]}"
            organ = parts[3]
            region_type = parts[4]  # enhancers or promoters
            
            data.append([species_direction, organ, region_type, output_dir])
        
        print(tabulate(data, headers=headers, tablefmt="grid"))
    
    return phase1_success and phase2_success 
