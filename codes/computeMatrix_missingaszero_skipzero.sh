#!/bin/sh


cd $dir

computeMatrix reference-point -o $output.matrix.gz --outFileNameMatrix $output.matrix.tab --outFileSortedRegions $output.regions.bed -S $scores -R $region -b $bef -a $bef --samplesLabel $samplesLabel --referencePoint center -bs $binSize --missingDataAsZero --skipZeros --blackListFileName $blacklist -p 5
