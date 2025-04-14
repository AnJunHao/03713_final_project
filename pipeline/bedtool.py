from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
from typing import NamedTuple
from pipeline.halper import HalperOutput

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
            assert halper_file.exists(), f"HALPER file {halper_file} does not exist"
        
        # Create output directory if it doesn't exist
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)

def load_bedtool_config(config_path: Path, halper_output: HalperOutput) -> BedtoolConfig:
    """
    Load bedtool configuration from a YAML file and a HalperOutput object.
    
    Args:
        config_path: Path to the config YAML file
        halper_output: A HalperOutput object containing the mapping file paths
        
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
        species_1_organ_1_to_species_2=halper_output.species_1_organ_1_to_species_2,
        species_1_organ_2_to_species_2=halper_output.species_1_organ_2_to_species_2,
        species_2_organ_1_to_species_1=halper_output.species_2_organ_1_to_species_1,
        species_2_organ_2_to_species_1=halper_output.species_2_organ_2_to_species_1,
    )

script_template = """#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --error={error_log}
#SBATCH --output={output_log}

module load bedtools

echo "Starting region comparison for {species_1} and {species_2}"

mkdir -p {output_dir}

# Path variables
SPECIES_1_ORGAN_1="{species_1_organ_1_peak_file}"
SPECIES_1_ORGAN_2="{species_1_organ_2_peak_file}"
SPECIES_2_ORGAN_1="{species_2_organ_1_peak_file}"
SPECIES_2_ORGAN_2="{species_2_organ_2_peak_file}"

# HALPER mappings: species_1 to species_2
SPECIES_1_ORGAN_1_TO_SPECIES_2="{species_1_organ_1_to_species_2}"
SPECIES_1_ORGAN_2_TO_SPECIES_2="{species_1_organ_2_to_species_2}"

# HALPER mappings: species_2 to species_1
SPECIES_2_ORGAN_1_TO_SPECIES_1="{species_2_organ_1_to_species_1}"
SPECIES_2_ORGAN_2_TO_SPECIES_1="{species_2_organ_2_to_species_1}"

echo "Comparing within-species (cross-tissue) overlaps..."

bedtools intersect -a "$SPECIES_1_ORGAN_1" -b "$SPECIES_1_ORGAN_2" -u > {output_dir}/{species_1}_shared_peaks.bed
bedtools intersect -a "$SPECIES_1_ORGAN_1" -b "$SPECIES_1_ORGAN_2" -v > {output_dir}/{species_1}_{organ_1}_specific.bed
bedtools intersect -a "$SPECIES_1_ORGAN_2" -b "$SPECIES_1_ORGAN_1" -v > {output_dir}/{species_1}_{organ_2}_specific.bed

bedtools intersect -a "$SPECIES_2_ORGAN_1" -b "$SPECIES_2_ORGAN_2" -u > {output_dir}/{species_2}_shared_peaks.bed
bedtools intersect -a "$SPECIES_2_ORGAN_1" -b "$SPECIES_2_ORGAN_2" -v > {output_dir}/{species_2}_{organ_1}_specific.bed
bedtools intersect -a "$SPECIES_2_ORGAN_2" -b "$SPECIES_2_ORGAN_1" -v > {output_dir}/{species_2}_{organ_2}_specific.bed

echo "Comparing cross-species conservation (species_1 to species_2)..."

bedtools intersect -a "$SPECIES_1_ORGAN_1_TO_SPECIES_2" -b "$SPECIES_2_ORGAN_1" -u > {output_dir}/{species_1}_{organ_1}_to_{species_2}_open.bed
bedtools intersect -a "$SPECIES_1_ORGAN_1_TO_SPECIES_2" -b "$SPECIES_2_ORGAN_1" -v > {output_dir}/{species_1}_{organ_1}_to_{species_2}_closed.bed

bedtools intersect -a "$SPECIES_1_ORGAN_2_TO_SPECIES_2" -b "$SPECIES_2_ORGAN_2" -u > {output_dir}/{species_1}_{organ_2}_to_{species_2}_open.bed
bedtools intersect -a "$SPECIES_1_ORGAN_2_TO_SPECIES_2" -b "$SPECIES_2_ORGAN_2" -v > {output_dir}/{species_1}_{organ_2}_to_{species_2}_closed.bed

echo "Comparing cross-species conservation (species_2 to species_1)..."

bedtools intersect -a "$SPECIES_2_ORGAN_1_TO_SPECIES_1" -b "$SPECIES_1_ORGAN_1" -u > {output_dir}/{species_2}_{organ_1}_to_{species_1}_open.bed
bedtools intersect -a "$SPECIES_2_ORGAN_1_TO_SPECIES_1" -b "$SPECIES_1_ORGAN_1" -v > {output_dir}/{species_2}_{organ_1}_to_{species_1}_closed.bed

bedtools intersect -a "$SPECIES_2_ORGAN_2_TO_SPECIES_1" -b "$SPECIES_1_ORGAN_2" -u > {output_dir}/{species_2}_{organ_2}_to_{species_1}_open.bed
bedtools intersect -a "$SPECIES_2_ORGAN_2_TO_SPECIES_1" -b "$SPECIES_1_ORGAN_2" -v > {output_dir}/{species_2}_{organ_2}_to_{species_1}_closed.bed

echo "Summary:"
for f in {output_dir}/*.bed; do
    echo "$(basename $f): $(wc -l < $f) peaks"
done

echo "Job finished"
"""

class GeneratedScriptOutput(NamedTuple):
    script: Path
    output_log: Path
    error_log: Path

def generate_script(config: BedtoolConfig) -> GeneratedScriptOutput:
    """
    Generate a script to compare peak regions between species and organs.

    Args:
        config: A BedtoolConfig object.
        temp_dir: Directory to store the generated script.

    Returns:
        A GeneratedScriptOutput containing paths to the script and log files.
    """
    temp_dir = config.temp_dir
    
    # Create script path and log paths
    script_name = f"compare_{config.species_1}_{config.species_2}_regions"
    script_path = temp_dir / f"{script_name}.job"
    error_log = config.output_dir / f"{script_name}.err.txt"
    output_log = config.output_dir / f"{script_name}.out.txt"
    
    # Delete old logs if they exist
    if output_log.exists():
        output_log.unlink()
    if error_log.exists():
        error_log.unlink()
    
    # Write the script
    with open(script_path, "w") as f:
        f.write(script_template.format(
            species_1=config.species_1,
            species_2=config.species_2,
            organ_1=config.organ_1,
            organ_2=config.organ_2,
            output_dir=config.output_dir,
            species_1_organ_1_peak_file=config.species_1_organ_1_peak_file,
            species_1_organ_2_peak_file=config.species_1_organ_2_peak_file,
            species_2_organ_1_peak_file=config.species_2_organ_1_peak_file,
            species_2_organ_2_peak_file=config.species_2_organ_2_peak_file,
            species_1_organ_1_to_species_2=config.species_1_organ_1_to_species_2,
            species_1_organ_2_to_species_2=config.species_1_organ_2_to_species_2,
            species_2_organ_1_to_species_1=config.species_2_organ_1_to_species_1,
            species_2_organ_2_to_species_1=config.species_2_organ_2_to_species_1,
            error_log=error_log,
            output_log=output_log
        ))
    
    return GeneratedScriptOutput(script_path, output_log, error_log)

def run_bedtool_pipeline(config_path: Path, halper_output) -> None:
    """
    Run the bedtools comparison pipeline.

    Args:
        config_path: Path to the configuration file.
        halper_output: HalperOutput object containing mapping files.
        temp_dir: Directory to store the generated script.
    """
    config = load_bedtool_config(config_path, halper_output)
    script_output = generate_script(config)
    script_path = script_output.script
    
    print(f"Submitting job: {script_path}")
    result = subprocess.run(["sbatch", str(script_path)], check=True, capture_output=True, text=True)
    
    if result.stdout:
        print(f"Job submission output: {result.stdout}")
    if result.stderr:
        print(f"Job submission error: {result.stderr}")
    
    print(f"Job submitted. Check logs at {script_output.output_log} and {script_output.error_log}")
