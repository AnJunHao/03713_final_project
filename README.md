# 03713 project



## 1. Install Required Tools

Follow the instructions in the official HAL liftover postprocessing guide to install HAL, HDF5, sonLib

🔗 [HAL installation instructions](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)



------

### 2. Mapping of chromatin regions cross species with halLiftover and HALPER (Human to mouse)


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

### Run `halLiftover` via `sbatch` Script

```bash
interact -n 4 -t 01:00:00 --mem=32gb
```

####  Pancreas Job Script (`human_to_mouse_pancreas.job`)

```bash
#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4

source ~/.bashrc
conda activate hal

bash /ocean/projects/bio230007p/yma6/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/yma6/data/idr.optimal_peak_pancreas.narrowPeak \
  -o /ocean/projects/bio230007p/yma6/halper_result \
  -s Human \
  -t Mouse \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
```

Submit the job:

```bash
sbatch human_to_mouse_pancreas.job
```

------

####  Liver Job Script (`human_to_mouse_liver.job`)

```bash
#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4

source ~/.bashrc
conda activate hal

bash /ocean/projects/bio230007p/yma6/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/yma6/data/idr.optimal_peak_liver.narrowPeak \
  -o /ocean/projects/bio230007p/yma6/halper_result \
  -s Human \
  -t Mouse \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
```

Submit the job:

```bash
sbatch human_to_mouse_liver.job
```


------



## 3.  Cross-Species and Cross-Tissue Comparison of Open Chromatin Regions

```
module load bedtools/2.30.0  
touch compare_open_regions.sh
```



```bash
#!/bin/bash

mkdir -p results

HUMAN_PANCREAS="data/idr.optimal_peak_pancreas_human.narrowPeak"
HUMAN_LIVER="data/idr.optimal_peak_liver_human.narrowPeak"
MOUSE_PANCREAS="data/idr.optimal_peak_pancreas_mouse.narrowPeak"
MOUSE_LIVER="data/idr.optimal_peak_liver_mouse.narrowPeak"

# HALPER output
HP_TO_MOUSE="halper_result/idr.optimal_peak_pancreas.HumanToMouse.HALPER.narrowPeak"
HL_TO_MOUSE="halper_result/idr.optimal_peak_liver.HumanToMouse.HALPER.narrowPeak"

echo "Comparing within-species (cross-tissue) overlaps..."

bedtools intersect -a "$HUMAN_PANCREAS" -b "$HUMAN_LIVER" -u > results/human_shared_peaks.bed
bedtools intersect -a "$HUMAN_PANCREAS" -b "$HUMAN_LIVER" -v > results/human_pancreas_specific.bed
bedtools intersect -a "$HUMAN_LIVER" -b "$HUMAN_PANCREAS" -v > results/human_liver_specific.bed

bedtools intersect -a "$MOUSE_PANCREAS" -b "$MOUSE_LIVER" -u > results/mouse_shared_peaks.bed
bedtools intersect -a "$MOUSE_PANCREAS" -b "$MOUSE_LIVER" -v > results/mouse_pancreas_specific.bed
bedtools intersect -a "$MOUSE_LIVER" -b "$MOUSE_PANCREAS" -v > results/mouse_liver_specific.bed

echo "Comparing cross-species conservation..."

bedtools intersect -a "$HP_TO_MOUSE" -b "$MOUSE_PANCREAS" -u > results/human_pancreas_to_mouse_open.bed
bedtools intersect -a "$HP_TO_MOUSE" -b "$MOUSE_PANCREAS" -v > results/human_pancreas_to_mouse_closed.bed

bedtools intersect -a "$HL_TO_MOUSE" -b "$MOUSE_LIVER" -u > results/human_liver_to_mouse_open.bed
bedtools intersect -a "$HL_TO_MOUSE" -b "$MOUSE_LIVER" -v > results/human_liver_to_mouse_closed.bed

echo "Summary:"
for f in results/*.bed; do
    echo "$(basename $f): $(wc -l < $f) peaks"
done

> /ocean/projects/bio230007p/yma6/halper_result/idr.optimal_peak_pancreas.HumanToMouse.HALPER.narrowPeak
```

```
Summary:
human_liver_specific.bed: 60831 peaks
human_liver_to_mouse_closed.bed: 41778 peaks
human_liver_to_mouse_open.bed: 37002 peaks
human_pancreas_specific.bed: 62075 peaks
human_pancreas_to_mouse_closed.bed: 54012 peaks
human_pancreas_to_mouse_open.bed: 43066 peaks
human_shared_peaks.bed: 122707 peaks
mouse_liver_specific.bed: 67512 peaks
mouse_pancreas_specific.bed: 25967 peaks
mouse_shared_peaks.bed: 55612 peaks
```













