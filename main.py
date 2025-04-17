#! /usr/bin/env python3

import pipeline
from pathlib import Path
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run HALPER and bedtools comparison pipeline")
    parser.add_argument(
        "--config", type=Path, default=Path("config.yaml"),
        help="Path to the config .yaml file (default: config.yaml)")
    args = parser.parse_args()
    
    print("Running HALPER pipeline...")
    success = pipeline.run_halper_pipeline(args.config, do_not_submit=True)
    print("HALPER pipeline complete!")
    # print(f"HALPER output files: {halper_output}")
    
    # print("\nRunning bedtools comparison pipeline...")
    # pipeline.run_bedtool_pipeline(args.config, halper_output)
    # print("Bedtools comparison pipeline submitted!")