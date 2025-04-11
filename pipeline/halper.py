from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml

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
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)

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
        output_dir=Path(config["output_dir"])
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

echo "Done"
"""

def generate_script(config: HalperConfig) -> Path:
    """
    Generate a script to map the peaks of the two species to each other.

    Args:
        config: A HalperConfig object.

    Returns:
        The path to the master script.
    """

    # We will do four mappings:
    # 1. species_1_organ_1_peak_file: species_1 -> species_2
    # 2. species_1_organ_2_peak_file: species_1 -> species_2
    # 3. species_2_organ_1_peak_file: species_2 -> species_1
    # 4. species_2_organ_2_peak_file: species_2 -> species_1
    
    script_paths = []
    
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
    
    # Create master script to submit all jobs
    master_script = config.temp_dir / "submit_all_jobs.sh"
    with open(master_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Submitting all jobs...'\n")
        for script in script_paths:
            f.write(f"sbatch {script}\n")
        f.write("echo 'All jobs submitted'\n")
    master_script.chmod(0o755)  # Make the master script executable
    
    return master_script

def run_halper_pipeline(config_path: Path) -> None:
    """
    Run the HALPER pipeline.

    Args:
        config: A HalperConfig object.
    
    Returns:
        None
    """
    config = load_halper_config(config_path)
    master_script = generate_script(config)
    result = subprocess.run(["bash", str(master_script)], check=True, capture_output=True, text=True)
    print(f"{result.stdout}")
    print(f"{result.stderr}")