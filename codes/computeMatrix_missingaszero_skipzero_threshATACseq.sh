#!/bin/sh

# a maximum threshold is used to remove effect of 1% outliers in the ATAC-seq signal

cd $dir

computeMatrix reference-point -o $output.matrix.gz --outFileNameMatrix $output.matrix.tab --outFileSortedRegions $output.regions.bed -S $scores -R $region -b $bef -a $bef --samplesLabel $samplesLabel --referencePoint center -bs $binSize --missingDataAsZero --skipZeros --blackListFileName $blacklist -p 5 ----maxThreshold 1000
