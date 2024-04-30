cd $dir

plotHeatmap -m $output.matrix.gz -out $output.heatmap.pdf --refPointLabel $refPointLabel --heatmapHeight 14 --colorMap $colorMap --yAxisLabel $yAxisLabel --xAxisLabel $xAxisLabel --regionsLabel $regionsLabel --zMin $zMin --zMax $zMax --yMin $yMin --yMax $yMax --outFileSortedRegions $output.sorted.regions.bed --sortRegions keep --averageTypeSummaryPlot mean --interpolationMethod bilinear
