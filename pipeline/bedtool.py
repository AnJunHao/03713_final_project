from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
from typing import NamedTuple
from pipeline.monitor import monitor_jobs

@dataclass
class BedtoolConfig:
    species_1: str
    species_2: str
    organ_1: str
    organ_2: str
    output_dir: Path
    temp_dir: Path
    
    # Peak files
    species_1_organ_1_peak_file: Path
    species_1_organ_2_peak_file: Path
    species_2_organ_1_peak_file: Path
    species_2_organ_2_peak_file: Path
    
    # HALPER output files - all four directions
    species_1_organ_1_to_species_2: Path
    species_1_organ_2_to_species_2: Path
    species_2_organ_1_to_species_1: Path
    species_2_organ_2_to_species_1: Path
    
    def __post_init__(self):
        # Check if peak files exist
        for peak_file in [self.species_1_organ_1_peak_file,
                          self.species_1_organ_2_peak_file,
                          self.species_2_organ_1_peak_file,
                          self.species_2_organ_2_peak_file]:
            assert peak_file.exists(), f"Peak file {peak_file} does not exist"
        
        # Check if HALPER files exist
        for halper_file in [self.species_1_organ_1_to_species_2,
                           self.species_1_organ_2_to_species_2,
                           self.species_2_organ_1_to_species_1,
                           self.species_2_organ_2_to_species_1]:
            # If the narrowPeak file exists, if not, check the gzipped file
            if not halper_file.exists():
                zipped_halper_file = halper_file.with_suffix(".gz")
                assert zipped_halper_file.exists(), f"HALPER file {zipped_halper_file} or {halper_file} does not exist"
                # Unzip the file to the same directory
                subprocess.run(["gunzip", zipped_halper_file])
                print(f"Unzipped HALPER file {zipped_halper_file} to {halper_file}")
            assert halper_file.exists(), f"HALPER file {halper_file} does not exist"
        
        # Create output directory if it doesn't exist
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory {self.output_dir}")

def load_bedtool_config(config_path: Path) -> BedtoolConfig:
    """
    Load bedtool configuration from a YAML file and a HalperOutput object.
    
    Args:
        config_path: Path to the config YAML file
        
    Returns:
        A BedtoolConfig object
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    return BedtoolConfig(
        species_1=config["species_1"],
        species_2=config["species_2"],
        organ_1=config["organ_1"],
        organ_2=config["organ_2"],
        output_dir=Path(config["bedtool_output_dir"]),
        temp_dir=Path(config["temp_dir"]),
        species_1_organ_1_peak_file=Path(config["species_1_organ_1_peak_file"]),
        species_1_organ_2_peak_file=Path(config["species_1_organ_2_peak_file"]),
        species_2_organ_1_peak_file=Path(config["species_2_organ_1_peak_file"]),
        species_2_organ_2_peak_file=Path(config["species_2_organ_2_peak_file"]),
        # Use the HalperOutput object for mapping files
        species_1_organ_1_to_species_2=Path(config["species_1_organ_1_to_species_2"]),
        species_1_organ_2_to_species_2=Path(config["species_1_organ_2_to_species_2"]),
        species_2_organ_1_to_species_1=Path(config["species_2_organ_1_to_species_1"]),
        species_2_organ_2_to_species_1=Path(config["species_2_organ_2_to_species_1"]),
    )

def bedtool_preprocess(config_path: Path) -> None:
    """
    Preprocess the configuration file for the bedtools pipeline.
    """
    # 1. Create backup of original config file
    backup_path = Path(f"{config_path}.backup02")
    with open(config_path, "r") as src:
        with open(backup_path, "w") as dst:
            dst.write(src.read())
    
    # 2. Read the config file
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    # 3. Unzip the HALPER files if they are zipped
    for entry in ["species_1_organ_1_to_species_2",
                  "species_1_organ_2_to_species_2",
                  "species_2_organ_1_to_species_1",
                  "species_2_organ_2_to_species_1"]:
        halper_file = Path(config[entry])
        if not halper_file.exists():
            # Check if the unzipped file exists
            unzipped_halper_file = halper_file.with_suffix("")
            if unzipped_halper_file.exists():
                halper_file = unzipped_halper_file
            else:
                raise FileNotFoundError(f"HALPER file {halper_file} or {unzipped_halper_file} does not exist")
        if halper_file.suffix == ".gz":
            subprocess.run(["gunzip", halper_file])
            halper_file = halper_file.with_suffix("")
        # Update the config file with the unzipped file
        config[entry] = str(halper_file)
    
    # 4. Extract the first 3 columns of the HALPER files
    cleaned_dir = Path(config["bedtool_output_dir"]) / "cleaned"
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    for entry in ["species_1_organ_1_to_species_2",
                  "species_1_organ_2_to_species_2",
                  "species_2_organ_1_to_species_1",
                  "species_2_organ_2_to_species_1"]:
        halper_file = Path(config[entry])
        cleaned_file = cleaned_dir / halper_file.name
        with open(cleaned_file, "w") as outfile:
            subprocess.run(["cut", "-f1-3", str(halper_file)], stdout=outfile)
        # Update the config file with the cleaned file
        config[entry+"_cleaned"] = str(cleaned_file)
    
    # 5. Extract the first 3 columns of the peak files
    for entry in ["species_1_organ_1_peak_file",
                  "species_1_organ_2_peak_file",
                  "species_2_organ_1_peak_file",
                  "species_2_organ_2_peak_file"]:
        peak_file = Path(config[entry])
        cleaned_file = cleaned_dir / peak_file.name
        with open(cleaned_file, "w") as outfile:
            subprocess.run(["cut", "-f1-3", str(peak_file)], stdout=outfile)
        # Update the config file with the cleaned file
        config[entry+"_cleaned"] = str(cleaned_file)
    
    # 6. Write the updated config file
    with open(config_path, "w") as f:
        yaml.dump(config, f)

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
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
bedtools intersect -a {halper_file} -b {native_file} -u > $conserved_file

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

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
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
    
    # Define the four combinations we want to analyze
    combinations = [
        {
            "prefix": f"{config.species_1}_{config.organ_1}_to_{config.species_2}_{config.organ_1}",
            "halper_file": config.species_1_organ_1_to_species_2,
            "native_file": config.species_2_organ_1_peak_file
        },
        {
            "prefix": f"{config.species_1}_{config.organ_2}_to_{config.species_2}_{config.organ_2}",
            "halper_file": config.species_1_organ_2_to_species_2,
            "native_file": config.species_2_organ_2_peak_file
        },
        {
            "prefix": f"{config.species_2}_{config.organ_1}_to_{config.species_1}_{config.organ_1}",
            "halper_file": config.species_2_organ_1_to_species_1,
            "native_file": config.species_1_organ_1_peak_file
        },
        {
            "prefix": f"{config.species_2}_{config.organ_2}_to_{config.species_1}_{config.organ_2}",
            "halper_file": config.species_2_organ_2_to_species_1,
            "native_file": config.species_1_organ_2_peak_file
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
        error_logs=error_logs
    )

def run_bedtool_pipeline(config_path: Path) -> bool:
    """
    Run the bedtools comparison pipeline.

    Args:
        config_path: Path to the configuration file.
    """
    bedtool_preprocess(config_path)
    config = load_bedtool_config(config_path)
    script_output = generate_script(config)
    script_path = script_output.script

    # Clean old output err logs
    for log in script_output.output_logs + script_output.error_logs:
        if log.exists():
            log.unlink()
            print(f"Deleted old log file: {log}")
    
    print(f"Submitting jobs: {script_path}")
    result = subprocess.run(["bash", str(script_path)], check=True, capture_output=True, text=True)
    
    if result.stdout:
        print(f"Job submission output: {result.stdout}")
    if result.stderr:
        print(f"Job submission error: {result.stderr}")
        return False
    
    # Monitor the submitted jobs
    try:
        success = monitor_jobs(script_output.output_logs, script_output.error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Jobs will still run.")
        success = False
    
    return success

class CrossSpeciesOrthologOpenVsClosedResult(NamedTuple):
    lifted_total: int
    open_total: int
    closed_total: int

def cross_species_ortholog_open_vs_closed(
        halper_file: Path,
        native_file: Path,
        prefix: str,
        output_dir: Path
        ) -> CrossSpeciesOrthologOpenVsClosedResult:
    """
    Identify the open chromatin regions in a given species whose orthologs in the other species are open or closed.
    This function implements the same logic as the bash script template used for the slurm jobs.
    
    Args:
        halper_file: Path to the cleaned HALPER file
        native_file: Path to the cleaned native peak file
        prefix: Prefix for the output files
        output_dir: Path to the output directory
        
    Returns:
        A CrossSpeciesOrthologOpenVsClosedResult object
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Intersect for conserved (ortholog open)
    conserved_file = output_dir / f"{prefix}_conserved.bed"
    print(f"[STEP 1] Finding conserved peaks: {conserved_file}")
    with open(conserved_file, "w") as out:
        subprocess.run(["bedtools", "intersect", "-a", str(halper_file), "-b", str(native_file), "-u"], stdout=out)

    # Step 2: Intersect for non-conserved (ortholog closed)
    closed_file = output_dir / f"{prefix}_closed.bed"
    print(f"[STEP 2] Finding non-conserved peaks: {closed_file}")
    with open(closed_file, "w") as out:
        subprocess.run(["bedtools", "intersect", "-a", str(halper_file), "-b", str(native_file), "-v"], stdout=out)

    # Step 3: Count peaks
    print("[STEP 3] Peak counts:")
    def count_lines(file: Path) -> int:
        return int(subprocess.check_output(["wc", "-l", str(file)]).decode().split()[0])
    
    lifted_total = count_lines(halper_file)
    open_total = count_lines(conserved_file)
    closed_total = count_lines(closed_file)
    
    print(f"Total lifted peaks: {lifted_total}")
    print(f"Open peaks: {open_total}")
    print(f"Closed peaks: {closed_total}")

    return CrossSpeciesOrthologOpenVsClosedResult(
        lifted_total=lifted_total,
        open_total=open_total,
        closed_total=closed_total
    )