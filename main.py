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
        "--skip-halper", action="store_true", 
        help="Skip running the HALPER pipeline (step 1).")
    parser.add_argument(
        "--skip-infer-halper-config-update", action="store_true", 
        help="Skip updating the config with HALPER output (step 1).")
    parser.add_argument(
        "--skip-cross-species", action="store_true", 
        help="Skip running cross-species ortholog open vs closed pipeline (step 3)")
    parser.add_argument(
        "--skip-cross-tissues", action="store_true", 
        help="Skip running cross-tissues region shared vs species pipeline (step 4)")
    args = parser.parse_args()
    
    print("="*100)
    if not args.skip_halper:
        print("Step 1: Running HALPER pipeline...")
        success = pipeline.run_halper_pipeline(args.config)
        print("Step 1: HALPER pipeline complete!")
    else:
        print("Step 1: Skipped HALPER pipeline")
        if not args.skip_infer_halper_config_update:
            success = pipeline.run_halper_pipeline(args.config, do_not_submit=True)
        else:
            print("Step 1: Skipped updating the config with HALPER output")

    print("="*100)
    print("Step 2: Preprocessing files for bedtools pipeline...")
    pipeline.bedtool_preprocess(args.config)
    print("Step 2: Bedtools preprocess complete!")
    
    print("="*100)
    if not args.skip_cross_species:
        print("Step 3: Running cross-species ortholog open vs closed pipeline...")
        success = pipeline.run_cross_species_ortholog_open_vs_closed_pipeline(args.config)
        print("Step 3: Cross-species ortholog open vs closed pipeline complete!")
    else:
        print("Step 3: Skipped cross-species ortholog open vs closed pipeline")

    print("="*100)
    if not args.skip_cross_tissues:
        print("Step 4: Running cross-tissues region shared vs species pipeline...")
        success = pipeline.run_cross_tissues_region_shared_vs_specific_pipeline(args.config)
        print("Step 4: Cross-tissues region shared vs species pipeline complete!")
    else:
        print("Step 4: Skipped cross-tissues region shared vs species pipeline")
