#!/bin/bash

# 输出文件夹
mkdir -p results

# 输入文件路径（你可以根据自己的路径进行修改）
HUMAN_PANCREAS="data/idr.optimal_peak_pancreas_human.narrowPeak"
HUMAN_LIVER="data/idr.optimal_peak_liver_human.narrowPeak"
MOUSE_PANCREAS="data/idr.optimal_peak_pancreas_mouse.narrowPeak"
MOUSE_LIVER="data/idr.optimal_peak_liver_mouse.narrowPeak"

# HALPER 输出文件（已经从 .gz 解压过）
HP_TO_MOUSE="halper_result/idr.optimal_peak_pancreas.HumanToMouse.HALPER.narrowPeak"
HL_TO_MOUSE="halper_result/idr.optimal_peak_liver.HumanToMouse.HALPER.narrowPeak"

# 1. 同物种跨组织 overlap
echo "Comparing within-species (cross-tissue) overlaps..."

bedtools intersect -a "$HUMAN_PANCREAS" -b "$HUMAN_LIVER" -u > results/human_shared_peaks.bed
bedtools intersect -a "$HUMAN_PANCREAS" -b "$HUMAN_LIVER" -v > results/human_pancreas_specific.bed
bedtools intersect -a "$HUMAN_LIVER" -b "$HUMAN_PANCREAS" -v > results/human_liver_specific.bed

bedtools intersect -a "$MOUSE_PANCREAS" -b "$MOUSE_LIVER" -u > results/mouse_shared_peaks.bed
bedtools intersect -a "$MOUSE_PANCREAS" -b "$MOUSE_LIVER" -v > results/mouse_pancreas_specific.bed
bedtools intersect -a "$MOUSE_LIVER" -b "$MOUSE_PANCREAS" -v > results/mouse_liver_specific.bed

# 3. 统计区域数量
echo "Summary:"
for f in results/*.bed; do
    echo "$(basename $f): $(wc -l < $f) peaks"
done
