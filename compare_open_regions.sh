#!/bin/bash

mkdir -p results

HUMAN_PANCREAS="data/idr.optimal_peak_pancreas_human.narrowPeak"
HUMAN_LIVER="data/idr.optimal_peak_liver_human.narrowPeak"
MOUSE_PANCREAS="data/idr.optimal_peak_pancreas_mouse.narrowPeak"
MOUSE_LIVER="data/idr.optimal_peak_liver_mouse.narrowPeak"

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
