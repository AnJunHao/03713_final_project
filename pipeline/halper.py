from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
from typing import NamedTuple, List
import time
from datetime import datetime

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
        if output_log.exists():
            output_log.unlink()
        if error_log.exists():
            error_log.unlink()
        
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

def monitor_jobs(output_logs: List[Path], error_logs: List[Path]) -> bool:
    """
    Monitor job status by checking output and error logs.
    
    Args:
        output_logs: List of paths to output log files
        error_logs: List of paths to error log files
        
    Returns:
        True when all jobs have completed successfully
    """
    job_statuses = {str(log): "QUEUED" for log in output_logs}
    all_complete = False
    start_time = datetime.now()
    
    # Print header once
    print("Monitoring job status...")
    
    # Initialize the status lines but don't print them yet
    status_lines = []
    for _ in range(len(output_logs)):
        status_lines.append("")
    
    # Add an extra line for elapsed time
    status_lines.append("")
    
    while not all_complete:
        all_complete = True
        
        # Clear previous status lines if any exist
        if status_lines[0]:  # If we've already printed status lines
            for _ in range(len(output_logs) + 1):  # +1 for elapsed time
                print("\033[A\033[K", end="")
        
        # Check and update each job's status
        for i, (out_log, err_log) in enumerate(zip(output_logs, error_logs)):
            job_name = out_log.name.replace(".out.txt", "")
            status_icon = "📋"
            status_text = "QUEUED"
            detail_text = ""
            
            # Check if logs exist
            if not out_log.exists() and not err_log.exists():
                job_statuses[str(out_log)] = "QUEUED"
                all_complete = False
            
            # Check if error log has content
            elif err_log.exists() and err_log.stat().st_size > 0:
                with open(err_log, 'r') as f:
                    err_content = f.read().strip()
                    if err_content:
                        job_statuses[str(out_log)] = "ERROR"
                        status_icon = "❌"
                        status_text = "ERROR"
                        # Get first line or first 50 chars of error
                        err_lines = err_content.splitlines()
                        detail_text = err_lines[-1][:50] + ("..." if len(err_lines[-1]) > 50 else "")
            
            # Check if output log exists and job is complete
            elif out_log.exists():
                if job_statuses[str(out_log)] == "QUEUED":
                    job_statuses[str(out_log)] = "RUNNING"
                
                # Check if job is complete by reading last line
                with open(out_log, 'r') as f:
                    content = f.read().strip()
                    if content:
                        lines = content.splitlines()
                        last_line = lines[-1] if lines else ""
                        
                        # Get last line of output for detail text
                        if lines:
                            detail_text = lines[-1][:50] + ("..." if len(lines[-1]) > 50 else "")
                        
                        if last_line == "Job finished":
                            job_statuses[str(out_log)] = "COMPLETED"
                            status_icon = "✅"
                            status_text = "COMPLETED"
                        else:
                            job_statuses[str(out_log)] = "RUNNING"
                            status_icon = "🏃"
                            status_text = "RUNNING"
                            all_complete = False
                    else:
                        job_statuses[str(out_log)] = "RUNNING"
                        status_icon = "🏃"
                        status_text = "RUNNING"
                        all_complete = False
            else:
                all_complete = False
            
            # Format and store the status line
            status_lines[i] = f"{status_icon} {job_name} {status_text}: {detail_text}"
            
            # Print current status
            print(status_lines[i])
        
        # Calculate and display elapsed time
        elapsed = datetime.now() - start_time
        hours, remainder = divmod(elapsed.total_seconds(), 3600)
        minutes, seconds = divmod(remainder, 60)
        elapsed_str = f"⏱️ Elapsed time: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"
        status_lines[-1] = elapsed_str
        print(elapsed_str)
        
        if not all_complete:
            time.sleep(1)  # Refresh every second
    
    # Count completed and error jobs
    completed_jobs = sum(1 for status in job_statuses.values() if status == "COMPLETED")
    error_jobs = sum(1 for status in job_statuses.values() if status == "ERROR")
    total_jobs = len(job_statuses)
    
    # Print final status with appropriate emoji
    if error_jobs == 0:
        print(f"\n🎉 All {total_jobs} jobs completed successfully in {elapsed_str[12:]}")
    elif completed_jobs == 0:
        print(f"\n❌ All {total_jobs} jobs failed")
    else:
        print(f"\n⚠️ {completed_jobs}/{total_jobs} jobs completed, {error_jobs} jobs failed")
    
    return error_jobs == 0

def run_halper_pipeline(config_path: Path, do_not_submit: bool = False) -> HalperOutput:
    """
    Run the HALPER pipeline.

    Args:
        config: A HalperConfig object.
    
    Returns:
        True if the pipeline was run successfully, False otherwise.
    """
    config = load_halper_config(config_path)
    script_output = generate_script(config)
    master_script = script_output.master_script
    output_logs = script_output.output_logs
    error_logs = script_output.error_logs
    halper_output = script_output.halper_output

    if do_not_submit:
        return halper_output

    result = subprocess.run(["bash", str(master_script)], check=True, capture_output=True, text=True)
    if result.stdout:
        print(f"{result.stdout}")
    if result.stderr:
        print(f"{result.stderr}")
        return halper_output
    
    # Monitor the submitted jobs
    try:
        monitor_jobs(output_logs, error_logs)
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Jobs will still run.")
    
    # Return the HALPER output
    return halper_output