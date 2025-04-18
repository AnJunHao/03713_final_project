# Cross-Species and Cross-Tissue Chromatin Accessibility Analysis Pipeline

This project provides a pipeline for comparing chromatin accessibility data (e.g., ATAC-seq peaks) across different species and tissues. It leverages HALPER for ortholog mapping and performs various downstream analyses to identify conserved, species-specific, and tissue-specific regulatory elements. Here is a 37-second demo video of the pipeline: 🎥 [Project Demo Video](https://youtu.be/aQjV1WrcCD0)

## 1. Setup

This pipeline is designed to be run on a HPC cluster and is tested on the [PSC Bridges-2 cluster](https://www.psc.edu/resources/bridges-2/). You can run the main script directly on the login node without any modification. The pipeline will automatically setup jobs on the compute nodes for each step.

### 1.1. Clone the Repository

```bash
git clone https://github.com/yitianTracy/03713_final_project.git
cd 03713_final_project
```

### 1.2. Python Environment

Ensure you have Python 3.9 or later installed.

```bash
python --version
```

Install the required Python packages:

```bash
pip install pyyaml tabulate
```

### 1.3. Install Required Tools

This pipeline relies on external bioinformatics tools. Please install them and ensure they are accessible in your system's PATH or provide the correct paths in the configuration file (`config.yaml`).

*   **HAL Tools:** Required for HALPER. Follow the official installation guide:
    *   🔗 [HAL installation instructions](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)
*   **BedTools:** Used for genomic interval operations. Installation instructions can be found on the [BedTools website](https://bedtools.readthedocs.io/en/latest/content/installation.html).

## 2. Configuration (`config.yaml`)

The pipeline's behavior is controlled by the `config.yaml` file. Before running, you need to configure it according to your environment and input data.

Key sections:

*   **Base Settings:** Define the species and organs being compared (`species_1`, `species_2`, `organ_1`, `organ_2`) and a temporary directory (`temp_dir`).
*   **HALPER Settings:**
    *   `halper_script`: Path to the `halper_map_peak_orthologs.sh` script from the [halLiftover-postprocessing](https://github.com/pfenninglab/halLiftover-postprocessing) tool.
    *   `hal_file`: Path to the HAL alignment file containing the genomes of interest.
    *   `halper_output_dir`: Directory where HALPER mapping results will be saved.
    *   `species_*_organ_*_peak_file`: Paths to the input peak files (e.g., narrowPeak format) for each species-organ combination.
*   **Enhancers vs Promoters Pipeline Settings:**
    *   `species_*_TSS_file`: Paths to the Transcription Start Site (TSS) annotation files (BED format) for each species. Used to distinguish promoter regions from enhancers.
*   **Output Directories:** Specify output directories for various analysis steps (`bedtool_preprocess_output_dir`, `cross_species_open_vs_closed_output_dir`, etc.). These directories will be created if they don't exist.

**Note:** The configuration file includes commented-out sections for paths that are *automatically inferred* by the pipeline (e.g., specific HALPER output files and conserved files). You typically do not need to uncomment or modify these.

Make sure all paths to input files, scripts, and tools are correct for your system.

## 3. Pipeline Execution

Run the entire pipeline using the main script:

```bash
python main.py --config config.yaml
```

Simply run the command in the **login node** and the pipeline will automatically setup jobs on the compute nodes for each step.

### Command-line Options

You can skip specific steps using command-line flags:

*   `--skip-step-1`: Skip running the HALPER pipeline (Step 1).
*   `--skip-infer-halper-output`: Skip updating the config with inferred HALPER output paths. This may be helpful if you want to manually configure the HALPER result files for subsequent steps.
*   `--skip-step-3`: Skip running the cross-species open vs. closed analysis (Step 3).
*   `--skip-infer-conserved-files`: Skip updating the config with inferred conserved file paths. This may be helpful if you want to manually configure the conserved files for subsequent steps.
*   `--skip-step-4`: Skip running the cross-tissues shared vs. specific analysis (Step 4).
*   `--skip-step-5`: Skip running the cross-tissues enhancer vs. promoter analysis (Step 5).
*   `--skip-step-6`: Skip running the cross-species enhancer vs. promoter analysis (Step 6).

**Note:** The output of step 1 is required for all subsequent steps. The output of step 3 is required for step 6. If you skip step 1 or step 3, you should satisfy one of the following conditions:

*   You have already run step 1 or step 3 before and have the output files saved in the corresponding directories without changing the names of the files. Our pipeline will automatically infer the paths to these output files and update the internal configuration for subsequent steps.
*   You have manually configured the required paths in the config file and use the `--skip-infer-halper-output` or `--skip-infer-conserved-files` flags.

Example: Run only Step 2 and Step 3 (assuming Step 1 outputs exist):

```bash
python main.py --config config.yaml --skip-step-1 --skip-step-4 --skip-step-5
```

## 4. Pipeline Walkthrough

The `main.py` script executes the following steps sequentially:

### Step 1: HALPER Ortholog Mapping (`pipeline/halper.py`)
Maps regulatory elements (peaks) from one species to another using the HALPER tool. This step:
- Takes narrowPeak files from each species-tissue combination
- Runs the HALPER ortholog mapping algorithm to find orthologous regions across species
- Generates mapped peak files in the target species' genome coordinates
- Creates the foundation for all downstream cross-species analyses

*Estimated runtime: ~3 hours*

### Step 2: BedTools Preprocessing (`pipeline/bedtool_preprocess.py`)
Prepares the data for the downstream analyses by:
- Extracting coordinate information (first 3 columns) from the peak files
- Unzipping HALPER output files if needed
- Creating clean BED files for use in subsequent steps
- Updating the configuration with paths to these cleaned files

*Estimated runtime: < 10 seconds*

### Step 3: Cross-Species Ortholog Comparison (Open vs. Closed) (`pipeline/cross_species_open_vs_closed.py`)
Identifies conserved and species-specific regulatory elements by:
- Comparing orthologous regions mapped from one species to corresponding native peaks in the target species
- Classifying orthologous elements as "open" (conserved) if they overlap with native peaks in the target species
- Classifying orthologous elements as "closed" (species-specific) if they don't overlap with native peaks
- Generating statistics on conservation rates between species for each tissue

*Estimated runtime (excluding queuing time): < 30 seconds*

### Step 4: Cross-Tissue Comparison (Shared vs. Specific) (`pipeline/cross_tissues_shared_vs_specific.py`)
Examines tissue specificity of regulatory elements within each species by:
- Comparing peak sets between different tissues of the same species
- Identifying "shared" peaks that appear in both tissues
- Identifying "tissue-specific" peaks unique to each tissue
- Generating statistics on the proportion of shared versus tissue-specific elements

*Estimated runtime (excluding queuing time): < 30 seconds*

### Step 5: Cross-Tissue Comparison (Enhancer vs. Promoter) (`pipeline/cross_tissues_enhancer_promoter.py`)
Classifies regulatory elements by genomic context and examines their tissue specificity by:
- Using TSS (Transcription Start Site) annotations to classify peaks as either promoters (≤5000bp from TSS) or enhancers (>5000bp from TSS)
- Comparing the distribution of enhancers and promoters between tissues of the same species
- Identifying tissue-shared versus tissue-specific enhancers and promoters
- Generating statistics on the proportion of enhancers versus promoters in each tissue

*Estimated runtime (excluding queuing time): < 30 seconds*

### Step 6: Cross-Species Comparison (Enhancer vs. Promoter) (`pipeline/cross_species_enhancer_promoter.py`)
Examines the evolutionary conservation of different types of regulatory elements by:
- Classifying conserved (open) peaks identified in Step 3 as enhancers or promoters based on distance to TSS
- Comparing enhancer and promoter conservation rates between species
- Identifying elements that maintain the same classification (enhancer/promoter) across species
- Generating statistics on conservation patterns specific to enhancers versus promoters

*Estimated runtime (excluding queuing time): < 30 seconds*

### Step 7: Motif Alignment

NOTE: this task is still a work in progress (the code/information associated with this task will be included in the final github repo/readme file)

### Additional Step: Functional Enrichment Analysis

We employed the GREAT (Genomic Regions Enrichment of Annotations Tool) to perform functional enrichment analysis of open chromatin regions. To run the enrichment analysis, users need to go to the [GREAT website](https://great.stanford.edu/) and upload corresponding BED files manually.

In our pipeline, we used the Mouse GRCm38 (UCSC mm10, December 2011) assembly as our reference genome. The output files from Step 4 were used as the input for the enrichment analysis to compare the enrichment of open chromatin regions that are conserved across tissues and specific to each tissue.

## 5. Outputs

Each pipeline step generates output files in the respective directories specified in `config.yaml`. Examine these directories after a successful run to find the results of each analysis (e.g., BED files, text tables). Error logs are also saved in the respective directories, with names ending in `*.err.txt`.

## 6. Additional Resources

*   [halLiftover-postprocessing Repository](https://github.com/pfenninglab/halLiftover-postprocessing)
*   [BedTools Documentation](https://bedtools.readthedocs.io/)
*   [HAL Tools](https://github.com/ComparativeGenomicsToolkit/hal)

## 7. Contributors

*   [Claude](https://github.com/AnJunHao)
*   [Isabella](https://github.com/iasalasallende)
*   [Shruthi](https://github.com/shruthirajaraman)
*   [Yitian](https://github.com/yitianTracy)