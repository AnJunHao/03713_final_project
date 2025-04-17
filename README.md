# 03713 project

## 1. Setup

### Clone the repository

```bash
git clone https://github.com/yitianTracy/03713_final_project.git
```

### Python Environment

Make sure your python version is at least 3.9.

```bash
python --version
```

Install the `pyyaml` package.

```bash
pip install pyyaml
```

### Install Required Tools

Follow the instructions in the official HAL liftover postprocessing guide to install HAL, HDF5, sonLib

🔗 [HAL installation instructions](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)

------

## 2. Mapping of chromatin regions cross species with halLiftover and HALPER (Human to mouse)

 🔗 [halLiftover-postprocessing](https://github.com/pfenninglab/halLiftover-postprocessing/tree/master)

### Select Human Peak Files

Chose the **optimal peak set** for human pancreas `/ocean/projects/bio230007p/ikaplow/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz` due to ~3-fold difference in sequencing depth (138M vs. 426M reads) between replicates.
 The **optimal set** is preferred over the conservative one for more comprehensive and reliable peak representation.

Chose the **optimal peak set** for human liver
 `/ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz` due to a moderate imbalance in sequencing depth and data quality between the two replicates. Replicate 1 shows significantly higher mitochondrial content (34.8% vs. 14.2%), a higher duplication rate (37.3% vs. 19.5%), anda ~2-fold difference in sequencing depth (95M vs. 201M reads). 

------

### Prepare Peak Files

```bash
# Pancreas (Human)
gunzip -c /ocean/projects/bio230007p/ikaplow/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
> idr.optimal_peak_pancreas.narrowPeak

# Liver (Human)
gunzip -c /ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz \
> idr.optimal_peak_liver.narrowPeak
```

------

### Configure the `config.yaml` file

```yaml
# Base Settings
species_1: "Human"
species_2: "Mouse"
organ_1: "Pancreas"
organ_2: "Liver"
temp_dir: "/ocean/projects/bio230007p/can1/temp"

# HALPER Settings
halper_script: "/ocean/projects/bio230007p/can1/halLiftover-postprocessing/halper_map_peak_orthologs.sh"
hal_file: "/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal"
halper_output_dir: "/ocean/projects/bio230007p/can1/output/halper"
species_1_organ_1_peak_file: "/ocean/projects/bio230007p/can1/s2_v2/human_pancreas.idr.optimal_peak.narrowPeak"
species_1_organ_2_peak_file: "/ocean/projects/bio230007p/can1/s2_v2/human_liver.idr.optimal_peak.narrowPeak"
species_2_organ_1_peak_file: "/ocean/projects/bio230007p/can1/s2_v2/mouse_pancreas.idr.optimal_peak.narrowPeak"
species_2_organ_2_peak_file: "/ocean/projects/bio230007p/can1/s2_v2/mouse_liver.idr.optimal_peak.narrowPeak"

# Bedtool Settings
bedtool_output_dir: "/ocean/projects/bio230007p/can1/output/bedtool"
```

Adjust the configurations to your own needs. Make sure the paths to the files / executables are correct. The `temp_dir` and `output_dir` will be created if they do not exist.

------

### Run the pipeline

```bash
python main.py --config config.yaml
```

The pipeline will run the HALPER pipeline and the bedtools pipeline in one go.

------

## 3.  Cross-Species and Cross-Tissue Comparison of Open Chromatin Regions

TODO