# Bin Comparisons
This file will contain documentation for each script used in the creation and analysis of copy number matrices and trees of varying bin sizes. Each script is designed for use on Linux OS, and use on other systems may not work.
## Dependencies
Each script contained assumes you have the following dependencies downloaded and properly setup:
* GetGiniSingle.r (as an environmental variable $GINISingle)
* GetGiniPaired.r (as an environmental variable $GINIPaired)
* SCOPESingle.r (as an environmental variable $SCOPESingle)
* SCOPEPaired.r (as an environmental variable $SCOPEPaired)
* SCOPE2CSV.py (as an environmental variable $SCOPE2CSV, located in PipelineScripts/Supplemental)
* S2M2.py (as an environmental variable $S2M2, located in PipelineScripts/Pipeline)
* [SeCNV](https://github.com/deepomicslab/SeCNV) (as an environmental variable $SECNV)
* [SCOPE](https://bioconductor.org/packages/release/bioc/html/SCOPE.html)
## Installation
To install,
```
git clone https://github.com/Carter-Payne/ResearchScripts/BinComparisonScripts.git
cd BinComparisonScripts
```
## Scripts
Below each scripts will be listed, along with a brief description:
* BinComparisons.py: Takes two folders of newick style unrooted trees and creates a .csv matrix containing each possible nRF tree comparison between the two folders.
* CNDistance.py: Returns the  distance of two different cell by location copy number matrices of equal bin length.
* FormatTSV.py: Takes the profiles.tsv file produced by CNAsim and formats it in a way to be properly used by CNDistance.py to measure the distance.
* GetGini(Single/Paired).r: produces the gini matrix for a given set of bam files at a given threshold for single/paired data.
* MultiCNA(Single/Paired).py: Creates multiple CNA matrices of varying bin sizes for a given set of bam files that are single/paired end.
* SCOPE(Single/Paired).r: Outputs a copy number and location matrix for a given bin size and gini threshold.
* TestKValues.py: Creates phylogenetic trees using SCOPE on different K values.
* PrintNormalCells.py: prints the inferred normal cells of a gini matrix by calculating a gini threshold.