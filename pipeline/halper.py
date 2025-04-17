from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
from typing import NamedTuple
from pipeline.monitor import monitor_jobs

@dataclass
class HalperConfig:
    species_1: str
    species_2: str
    organ_1: str
    organ_2: str
    temp_dir: Path
    halper_script: Path
    hal_file: Path
    species_1_organ_1_peak_file: Path
    species_1_organ_2_peak_file: Path
    species_2_organ_1_peak_file: Path
    species_2_organ_2_peak_file: Path
    output_dir: Path

    def __post_init__(self):
        # Check if peak files exist and are indeed narrowPeak files
        for peak_file in [self.species_1_organ_1_peak_file,
                          self.species_1_organ_2_peak_file,
                          self.species_2_organ_1_peak_file,
                          self.species_2_organ_2_peak_file]:
            assert peak_file.exists(), f"Peak file {peak_file} does not exist"
            assert peak_file.suffix == ".narrowPeak", f"Peak file {peak_file} is not a narrowPeak file"
        # Check if HAL file exists and is indeed a HAL file
        assert self.hal_file.exists(), f"HAL file {self.hal_file} does not exist"
        assert self.hal_file.suffix == ".hal", f"HAL file {self.hal_file} is not a HAL file"
        # Check if temp and output directories exist
        # Create them if they don't exist
        if not self.temp_dir.exists():
            self.temp_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created temp directory {self.temp_dir}")
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory {self.output_dir}")

def load_halper_config(config_path: Path) -> HalperConfig:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    return HalperConfig(
        species_1=config["species_1"],
        species_2=config["species_2"],
        organ_1=config["organ_1"],
        organ_2=config["organ_2"],
        temp_dir=Path(config["temp_dir"]),
        halper_script=Path(config["halper_script"]),
        hal_file=Path(config["hal_file"]),
        species_1_organ_1_peak_file=Path(config["species_1_organ_1_peak_file"]),
        species_1_organ_2_peak_file=Path(config["species_1_organ_2_peak_file"]),
        species_2_organ_1_peak_file=Path(config["species_2_organ_1_peak_file"]),
        species_2_organ_2_peak_file=Path(config["species_2_organ_2_peak_file"]),
        output_dir=Path(config["halper_output_dir"]) # NOTE: Different naming! Be careful!
    )

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load anaconda3
conda init
source ~/.bashrc

echo "mapping {source_species}: {source_organ} to {target_species}"

bash {halper_script} \
  -b {peak_file} \
  -o {output_dir} \
  -s {source_species} \
  -t {target_species} \
  -c {hal_file}

echo "Job finished"
"""

class HalperOutput(NamedTuple):
    species_1_organ_1_to_species_2: Path
    species_1_organ_2_to_species_2: Path
    species_2_organ_1_to_species_1: Path
    species_2_organ_2_to_species_1: Path

class GeneratedScriptOutput(NamedTuple):
    master_script: Path
    output_logs: list[Path]
    error_logs: list[Path]
    halper_output: HalperOutput

def generate_script(config: HalperConfig) -> GeneratedScriptOutput:
    """
    Generate a script to map the peaks of the two species to each other.

    Args:
        config: A HalperConfig object.

    Returns:
        A tuple containing:
        1. the path to the master script
        2. a list of paths to the output logs
        3. a list of paths to the error logs
    """

    # We will do four mappings:
    # 1. species_1_organ_1_peak_file: species_1 -> species_2
    # 2. species_1_organ_2_peak_file: species_1 -> species_2
    # 3. species_2_organ_1_peak_file: species_2 -> species_1
    # 4. species_2_organ_2_peak_file: species_2 -> species_1
    
    script_paths = []
    output_logs = []
    error_logs = []
    halper_output = []
    
    # Define mapping configurations
    mappings = [
        {
            "source_species": config.species_1, 
            "source_organ": config.organ_1, 
            "target_species": config.species_2,
            "peak_file": config.species_1_organ_1_peak_file
        },
        {
            "source_species": config.species_1, 
            "source_organ": config.organ_2, 
            "target_species": config.species_2,
            "peak_file": config.species_1_organ_2_peak_file
        },
        {
            "source_species": config.species_2, 
            "source_organ": config.organ_1, 
            "target_species": config.species_1,
            "peak_file": config.species_2_organ_1_peak_file
        },
        {
            "source_species": config.species_2, 
            "source_organ": config.organ_2, 
            "target_species": config.species_1,
            "peak_file": config.species_2_organ_2_peak_file
        }
    ]
    
    # Create scripts for each mapping
    for mapping in mappings:
        source_species = mapping["source_species"]
        source_organ = mapping["source_organ"]
        target_species = mapping["target_species"]
        peak_file = mapping["peak_file"]
        
        script = config.temp_dir / f"map_{source_species}_{source_organ}_to_{target_species}.job"
        error_log = config.output_dir / f"map_{source_species}_{source_organ}_to_{target_species}.err.txt"
        output_log = config.output_dir / f"map_{source_species}_{source_organ}_to_{target_species}.out.txt"

        with open(script, "w") as f:
            f.write(script_template.format(
                source_species=source_species,
                source_organ=source_organ,
                target_species=target_species,
                halper_script=config.halper_script,
                peak_file=peak_file,
                output_dir=config.output_dir,
                hal_file=config.hal_file,
                error_log=error_log,
                output_log=output_log
            ))
        script_paths.append(script)
        output_logs.append(output_log)
        error_logs.append(error_log)

        # Delete old output and error logs if they exist
        # if output_log.exists():
        #     output_log.unlink()
        # if error_log.exists():
        #     error_log.unlink()
        
        # Determine the HALPER output file name
        # Extract the base name from the peak file (without extension)
        peak_base_name = peak_file.stem
        # Format the mapping string (SourceToTarget format)
        mapping_string = f"{source_species}To{target_species}"
        # Construct the full output path with the fixed suffix
        halper_output_file = config.output_dir / f"{peak_base_name}.{mapping_string}.HALPER.narrowPeak.gz"
        halper_output.append(halper_output_file)

    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all jobs...'\n")
        for script in script_paths:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    master_script.chmod(0o755)  # Make the master script executable
    
    # Create the HalperOutput object with the expected output files
    halper_output_result = HalperOutput(
        species_1_organ_1_to_species_2=halper_output[0],
        species_1_organ_2_to_species_2=halper_output[1],
        species_2_organ_1_to_species_1=halper_output[2],
        species_2_organ_2_to_species_1=halper_output[3]
    )
    
    return GeneratedScriptOutput(master_script,
                                 output_logs,
                                 error_logs,
                                 halper_output_result)

def update_config(config_path: Path, halper_output: HalperOutput) -> None:
    """
    Update the config file with the HalperOutput entries.
    
    Args:
        config_path: Path to the config file
        halper_output: HalperOutput object containing output file paths
    """
    # Create backup of original config file
    backup_path = Path(f"{config_path}.backup01")
    with open(config_path, "r") as src:
        with open(backup_path, "w") as dst:
            dst.write(src.read())
    
    # Load existing config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    # Add HalperOutput entries to config
    config["species_1_organ_1_to_species_2"] = str(halper_output.species_1_organ_1_to_species_2)
    config["species_1_organ_2_to_species_2"] = str(halper_output.species_1_organ_2_to_species_2)
    config["species_2_organ_1_to_species_1"] = str(halper_output.species_2_organ_1_to_species_1)
    config["species_2_organ_2_to_species_1"] = str(halper_output.species_2_organ_2_to_species_1)
    
    # Write updated config back to file
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_flow_style=False)

def run_halper_pipeline(config_path: Path, do_not_submit: bool = False) -> bool:
    """
    Run the HALPER pipeline.

    Args:
        config: A HalperConfig object.
    
    Returns:
        True if all jobs completed successfully, False otherwise.
    """
    config = load_halper_config(config_path)
    script_output = generate_script(config)
    master_script = script_output.master_script
    output_logs = script_output.output_logs
    error_logs = script_output.error_logs
    halper_output = script_output.halper_output
    update_config(config_path, halper_output)

    if do_not_submit:
        print("HALPER pipeline is skipped!")
        return True
    
    # Delete old output and error logs if they exist
    old_log_count = 0
    for log in output_logs + error_logs:
        if log.exists():
            log.unlink()
            old_log_count += 1
    if old_log_count > 0:
        print(f"Deleted {old_log_count} old log files")

    # Submit the jobs
    result = subprocess.run(["bash", str(master_script)], check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return False
    
    # Monitor the submitted jobs
    try:
        success = monitor_jobs(output_logs, error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Jobs will still run.")
        success = False
    
    return success