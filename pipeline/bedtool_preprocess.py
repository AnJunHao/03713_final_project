from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml

@dataclass
class Config:
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
    
    # HALPER output files - all four directions (optional)
    species_1_organ_1_to_species_2: Path | None
    species_1_organ_2_to_species_2: Path | None
    species_2_organ_1_to_species_1: Path | None
    species_2_organ_2_to_species_1: Path | None
    
    # TSS files
    species_1_tss_file: Path
    species_2_tss_file: Path

    # Genome fasta files
    species_1_genome_fasta: Path
    species_2_genome_fasta: Path

    # Conserved files (optional)
    species_1_to_species_2_organ_1_conserved: Path | None
    species_1_to_species_2_organ_2_conserved: Path | None
    species_2_to_species_1_organ_1_conserved: Path | None
    species_2_to_species_1_organ_2_conserved: Path | None

    # Promoter/enhancer files (optional)
    species_1_organ_1_promoters: Path | None
    species_1_organ_1_enhancers: Path | None
    species_1_organ_2_promoters: Path | None
    species_1_organ_2_enhancers: Path | None
    species_2_organ_1_promoters: Path | None
    species_2_organ_1_enhancers: Path | None
    species_2_organ_2_promoters: Path | None
    species_2_organ_2_enhancers: Path | None

    # Mapped Promoter/enhancer files (optional)
    species_1_to_species_2_organ_1_conserved_promoters: Path | None
    species_1_to_species_2_organ_1_conserved_enhancers: Path | None
    species_1_to_species_2_organ_2_conserved_promoters: Path | None
    species_1_to_species_2_organ_2_conserved_enhancers: Path | None
    species_2_to_species_1_organ_1_conserved_promoters: Path | None
    species_2_to_species_1_organ_1_conserved_enhancers: Path | None
    species_2_to_species_1_organ_2_conserved_promoters: Path | None
    species_2_to_species_1_organ_2_conserved_enhancers: Path | None

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
            if halper_file is None:
                continue
            # If the narrowPeak file exists, if not, check the gzipped file
            if not halper_file.exists():
                zipped_halper_file = halper_file.with_suffix(".gz")
                assert zipped_halper_file.exists(), f"HALPER file {zipped_halper_file} or {halper_file} does not exist"
                # Unzip the file to the same directory
                subprocess.run(["gunzip", zipped_halper_file])
                print(f"Unzipped HALPER file {zipped_halper_file} to {halper_file}")
            assert halper_file.exists(), f"HALPER file {halper_file} does not exist"
        
        # Check if TSS files exist
        for tss_file in [self.species_1_tss_file,
                         self.species_2_tss_file]:
            assert tss_file.exists(), f"TSS file {tss_file} does not exist"
        
        # Create output directory if it doesn't exist
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory {self.output_dir}")

def load_bedtool_config(config_path: Path, output_dir_entry: str) -> Config:
    """
    Load bedtool configuration from a YAML file and a HalperOutput object.
    
    Args:
        config_path: Path to the config YAML file
        output_dir_entry: The config entry to use as output directory
        
    Returns:
        A BedtoolConfig object
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    return Config(
        species_1=config["species_1"],
        species_2=config["species_2"],
        organ_1=config["organ_1"],
        organ_2=config["organ_2"],
        output_dir=Path(config[output_dir_entry]),
        temp_dir=Path(config["temp_dir"]),
        species_1_organ_1_peak_file=Path(config["species_1_organ_1_peak_file_cleaned"]),
        species_1_organ_2_peak_file=Path(config["species_1_organ_2_peak_file_cleaned"]),
        species_2_organ_1_peak_file=Path(config["species_2_organ_1_peak_file_cleaned"]),
        species_2_organ_2_peak_file=Path(config["species_2_organ_2_peak_file_cleaned"]),
        # Use the HalperOutput object for mapping files
        species_1_organ_1_to_species_2=Path(config["species_1_organ_1_to_species_2_cleaned"])\
            if "species_1_organ_1_to_species_2_cleaned" in config else None,
        species_1_organ_2_to_species_2=Path(config["species_1_organ_2_to_species_2_cleaned"])\
            if "species_1_organ_2_to_species_2_cleaned" in config else None,
        species_2_organ_1_to_species_1=Path(config["species_2_organ_1_to_species_1_cleaned"])\
            if "species_2_organ_1_to_species_1_cleaned" in config else None,
        species_2_organ_2_to_species_1=Path(config["species_2_organ_2_to_species_1_cleaned"])\
            if "species_2_organ_2_to_species_1_cleaned" in config else None,
        species_1_tss_file=Path(config["species_1_TSS_file"]),
        species_2_tss_file=Path(config["species_2_TSS_file"]),
        species_1_genome_fasta=Path(config["species_1_genome_fasta"]),
        species_2_genome_fasta=Path(config["species_2_genome_fasta"]),
        # Conserved files (optional)
        species_1_to_species_2_organ_1_conserved=Path(config["species_1_to_species_2_organ_1_conserved"])\
            if "species_1_to_species_2_organ_1_conserved" in config else None,
        species_1_to_species_2_organ_2_conserved=Path(config["species_1_to_species_2_organ_2_conserved"])\
            if "species_1_to_species_2_organ_2_conserved" in config else None,
        species_2_to_species_1_organ_1_conserved=Path(config["species_2_to_species_1_organ_1_conserved"])\
            if "species_2_to_species_1_organ_1_conserved" in config else None,
        species_2_to_species_1_organ_2_conserved=Path(config["species_2_to_species_1_organ_2_conserved"])\
            if "species_2_to_species_1_organ_2_conserved" in config else None,
        # Promoter/enhancer files (optional)
        species_1_organ_1_promoters=Path(config["species_1_organ_1_promoters"])\
            if "species_1_organ_1_promoters" in config else None,
        species_1_organ_1_enhancers=Path(config["species_1_organ_1_enhancers"])\
            if "species_1_organ_1_enhancers" in config else None,
        species_1_organ_2_promoters=Path(config["species_1_organ_2_promoters"])\
            if "species_1_organ_2_promoters" in config else None,
        species_1_organ_2_enhancers=Path(config["species_1_organ_2_enhancers"])\
            if "species_1_organ_2_enhancers" in config else None,
        species_2_organ_1_promoters=Path(config["species_2_organ_1_promoters"])\
            if "species_2_organ_1_promoters" in config else None,   
        species_2_organ_1_enhancers=Path(config["species_2_organ_1_enhancers"])\
            if "species_2_organ_1_enhancers" in config else None,
        species_2_organ_2_promoters=Path(config["species_2_organ_2_promoters"])\
            if "species_2_organ_2_promoters" in config else None,
        species_2_organ_2_enhancers=Path(config["species_2_organ_2_enhancers"])\
            if "species_2_organ_2_enhancers" in config else None,
        # Mapped Promoter/enhancer files (optional)
        species_1_to_species_2_organ_1_conserved_promoters=Path(config["species_1_to_species_2_organ_1_conserved_promoters"])\
            if "species_1_to_species_2_organ_1_conserved_promoters" in config else None,
        species_1_to_species_2_organ_1_conserved_enhancers=Path(config["species_1_to_species_2_organ_1_conserved_enhancers"])\
            if "species_1_to_species_2_organ_1_conserved_enhancers" in config else None,
        species_1_to_species_2_organ_2_conserved_promoters=Path(config["species_1_to_species_2_organ_2_conserved_promoters"])\
            if "species_1_to_species_2_organ_2_conserved_promoters" in config else None,
        species_1_to_species_2_organ_2_conserved_enhancers=Path(config["species_1_to_species_2_organ_2_conserved_enhancers"])\
            if "species_1_to_species_2_organ_2_conserved_enhancers" in config else None,
        species_2_to_species_1_organ_1_conserved_promoters=Path(config["species_2_to_species_1_organ_1_conserved_promoters"])\
            if "species_2_to_species_1_organ_1_conserved_promoters" in config else None,
        species_2_to_species_1_organ_1_conserved_enhancers=Path(config["species_2_to_species_1_organ_1_conserved_enhancers"])\
            if "species_2_to_species_1_organ_1_conserved_enhancers" in config else None,
        species_2_to_species_1_organ_2_conserved_promoters=Path(config["species_2_to_species_1_organ_2_conserved_promoters"])\
            if "species_2_to_species_1_organ_2_conserved_promoters" in config else None,
        species_2_to_species_1_organ_2_conserved_enhancers=Path(config["species_2_to_species_1_organ_2_conserved_enhancers"])\
            if "species_2_to_species_1_organ_2_conserved_enhancers" in config else None,
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
        if entry not in config:
            continue
        halper_file = Path(config[entry])
        if not halper_file.exists():
            # Check if the unzipped file exists
            unzipped_halper_file = halper_file.with_suffix("")
            if unzipped_halper_file.exists():
                halper_file = unzipped_halper_file
            else:
                raise FileNotFoundError(f"HALPER file {halper_file} or {unzipped_halper_file} does not exist")
        if halper_file.suffix == ".gz":
            if not halper_file.with_suffix("").exists():
                subprocess.run(["gunzip", halper_file])
            halper_file = halper_file.with_suffix("")
        # Update the config file with the unzipped file
        config[entry] = str(halper_file)
    
    # 4. Extract the first 3 columns of the HALPER files
    cleaned_dir = Path(config["bedtool_preprocess_output_dir"])
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    for entry in ["species_1_organ_1_to_species_2",
                  "species_1_organ_2_to_species_2",
                  "species_2_organ_1_to_species_1",
                  "species_2_organ_2_to_species_1"]:
        if entry not in config:
            continue
        halper_file = Path(config[entry])
        cleaned_file = cleaned_dir / halper_file.name
        if not cleaned_file.exists():
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
        if not cleaned_file.exists():
            with open(cleaned_file, "w") as outfile:
                subprocess.run(["cut", "-f1-3", str(peak_file)], stdout=outfile)
        # Update the config file with the cleaned file
        config[entry+"_cleaned"] = str(cleaned_file)
    
    # 6. Write the updated config file
    with open(config_path, "w") as f:
        yaml.dump(config, f)