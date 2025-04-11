#! /usr/bin/env python3

import pipeline
from pathlib import Path
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run HALPER and bedtools comparison pipeline")
    parser.add_argument(
        "--config", type=Path, default=Path("config.yaml"),
        help="Path to the config .yaml file (default: config.yaml)")
    parser.add_argument(
        "--temp_dir", type=Path, default=Path("temp"),
        help="Path to the temporary directory (default: temp)")
    args = parser.parse_args()
    
    # Create temp directory if it doesn't exist
    if not args.temp_dir.exists():
        args.temp_dir.mkdir(parents=True, exist_ok=True)
    
    print("Running HALPER pipeline...")
    halper_output = pipeline.run_halper_pipeline(args.config)
    print("HALPER pipeline complete!")
    print(f"HALPER output files: {halper_output}")
    
    print("\nRunning bedtools comparison pipeline...")
    pipeline.run_bedtool_pipeline(args.config, halper_output, args.temp_dir)
    print("Bedtools comparison pipeline submitted!")