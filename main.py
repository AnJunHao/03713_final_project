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
        "--skip-step-1", action="store_true", 
        help="Skip running the HALPER pipeline (step 1).")
    parser.add_argument(
        "--skip-step-3", action="store_true", 
        help="Skip running cross-species ortholog open vs closed pipeline (step 3).")
    parser.add_argument(
        "--skip-step-4", action="store_true", 
        help="Skip running cross-tissues region shared vs species pipeline (step 4)")
    parser.add_argument(
        "--skip-step-5", action="store_true", 
        help="Skip running cross-tissues (within-species) enhancers vs promoters pipeline (step 5)")
    parser.add_argument(
        "--skip-step-6", action="store_true", 
        help="Skip running cross-species enhancers vs promoters pipeline (step 6)")
    parser.add_argument(
        "--skip-step-7", action="store_true", 
        help="Skip running MEME-chip analysis pipeline (step 7)")
    args = parser.parse_args()
    
    print("="*100)
    if not args.skip_step_1:
        print("Step 1: Running HALPER pipeline...")
        success = pipeline.run_halper_pipeline(args.config)
        print("Step 1: HALPER pipeline complete!")
    else:
        print("Step 1: Skipped HALPER pipeline")

    print("="*100)
    print("Step 2: Preprocessing files for bedtools pipeline...")
    pipeline.bedtool_preprocess(args.config)
    print("Step 2: Bedtools preprocess complete!")
    
    print("="*100)
    if not args.skip_step_3:
        print("Step 3: Running cross-species ortholog open vs closed pipeline...")
        success = pipeline.run_cross_species_open_vs_closed_pipeline(args.config)
        print("Step 3: Cross-species ortholog open vs closed pipeline complete!")
    else:
        print("Step 3: Skipped cross-species ortholog open vs closed pipeline")

    print("="*100)
    if not args.skip_step_4:
        print("Step 4: Running cross-tissues region shared vs species pipeline...")
        success = pipeline.run_cross_tissues_shared_vs_specific_pipeline(args.config)
        print("Step 4: Cross-tissues region shared vs species pipeline complete!")
    else:
        print("Step 4: Skipped cross-tissues region shared vs species pipeline")

    print("="*100)
    if not args.skip_step_5:
        print("Step 5: Running cross-tissues (within-species) enhancers vs promoters pipeline...")
        success = pipeline.run_cross_tissues_enhancer_promoter_pipeline(args.config)
        print("Step 5: Cross-tissues (within-species) enhancers vs promoters pipeline complete!")
    else:
        print("Step 5: Skipped cross-tissues (within-species) enhancers vs promoters pipeline")

    print("="*100)
    if not args.skip_step_6:
        print("Step 6: Running cross-species enhancers vs promoters pipeline...")
        success = pipeline.run_cross_species_enhancer_promoter_pipeline(args.config)
        print("Step 6: Cross-species enhancers vs promoters pipeline complete!")
    else:
        print("Step 6: Skipped cross-species enhancers vs promoters pipeline")

    print("="*100)
    if not args.skip_step_7:
        print("Step 7: Running MEME-chip analysis pipeline...")
        success = pipeline.run_meme_chip_analysis_pipeline(args.config)
        print("Step 7: MEME-chip analysis pipeline complete!")
    else:
        print("Step 7: Skipped MEME-chip analysis pipeline")
