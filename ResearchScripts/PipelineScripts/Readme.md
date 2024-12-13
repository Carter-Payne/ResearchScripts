# Sequencing Pipeline
This file will contain documentation for each script used in the pipeline and analysis of sequenced cancer data. Each script is designed for use on Linux OS, and use on other systems may not work.
## Dependencies
Each script contained assumes you have the following dependencies downloaded and added to path:
* [bwa](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools/samtools)
* [picard](https://github.com/broadinstitute/picard) (as an environmental variable $PICARD)
* [gatk](https://gatk.broadinstitute.org/hc/en-us)
* [bcftools](https://samtools.github.io/bcftools/)
* [bedtools](https://github.com/arq5x/bedtools2)
* [gnu parallel](https://www.gnu.org/software/parallel/)
* [SeCNV](https://github.com/deepomicslab/SeCNV) (as an environmental variable $SECNV)
## Installation
To install,
```
git clone https://github.com/Carter-Payne/ResearchScripts/PipelineScripts.git
cd PreprocessingScripts
```
## Pipeline.py
This is the main script, which utilizes FASTERQ.sh, Preprocess.py, SNV.py, CNA.py, and S2M2.py. This Pipeline is designed for single-cell single sample processing, so if you want to process multiple samples, it is advised to process each one individually.

Usage is the following:
```
python3 Pipeline.py [-h] -a Accession List -ao FASTQ folder -pd Prefetch folder -r Reference file [-i] [-t Threads] -rg ReadGroup data [ReadGroup data ...] [-pu] -temp Temporary folder [-rem] -o Bam output [-snv SNV folder] [-cna CNA folder] [-hg Genome type] [-bqsr Base Quality Recalibration [Base Quality Recalibration ...]] [-regions VCFRegions] [-mq Mapping Quality] [-pq Phred Quality] [-mind Minimum Read Depth] [-maxd Maximum Read Depth] -p File Prefix [-b Bin sizes] [-d] [-l Log file]
```
* -h, --help: display the help information
* -a, --ACCESSION: Path to the accession list of the sample(s) you want to process
* -ao, --ACCESSIONOUTPUT: Path to the output folder for the FASTQ files
* -pd, --PREFETCHDIR: Path to the output folder of the prefetch directory, set by the SRA toolkit
* -r, --REFERENCE: Path to the reference file
* -i, --INDEX: Set if you need to index the reference file
* -t, --THREADS: Number of threads to use (default: Auto)
* -rg, --READGROUP: Read group data for each run in the format of PL LB SAMN. If SAMN is not given, each run will be given its own SAMN
* -pg, --PLATFORMUNIT: select if you want each run to use the platform unit found in its FASTQ files(If it wasn't found, 'machine' will be used). Adding read groups will take longer if this is selected. If it is not, the Platform Unit will be the ID.
* -temp, --TEMP: path to temporary folder of the outputs of the pipeline
* -rem, --REMOVE: select if you wish to remove duplicates instead of marking them
* -o, --OUTPUT: Path to the finalized .bam files
* -snv, --SNVPATH: Path to the output folder for the finished .vcf file. If not given, .vcf file won't be created
* -cna, --CNAPATH: Path to the output folder for the .csv file. If not given, .csv file won't be created
* -h, --HG: The type of reference genome being used, either hg38 or hg19(default: hg38)
* -bqsr, --BASEQUALITY: Paths to each file to be used for BQSR. If none are given, BQSR will not be done.
* -regions, --REGIONS: Path to a text file containing which chromosomes to perform variant calling on(default: all chromosomes in reference file)
* -mq, --MQ: Minimum mapping quality for read filtering(default: 40)
* -pq, --PQ: Minimum phred quality for variant filtering (default: 40)
* -mind, --MINDEPTH: Minimum read depth for variant filtering(default: 10)
* -maxd, --MAXDEPTH: Maximum read depth for variant filtering(default: None)
* -p, --PREFIX: Prefix for the finished CNA and SNV files.
* -b, --BIN: Bin sizes to be used for CNA analysis(default: 500000)
* -d, --DEBUG: Select this option if you need to keep each intermitten file in the pipeline(only select this option if you have enough space to do so.)
* -l, --LOG: Path for the file storing pipeline running times.(default: current directory)

Example:
```
python3 Pipeline.py -a /path/to/AccessionList.txt -ao /path/to/AccessionOutput/folder -pd /path/to/SRA/Prefetch -r /path/to/referencefile.fa -i -t 10 -rg Illumina HiSeq2500 -temp /path/to/temp/folder -o /path/to/bam/output -snv /path/to/snv/folder -cna /path/to/cna/folder -bqsr /path/to/dbSNP.vcf.gz /path/to/AnotherdbSNP.vcf.gz -mq 40 -pq 40 -mind 20 -maxd 250 -p Example -b 500000 -l /path/to/log/folder -d
```
The accessory scripts to the main Pipeline.py script are the following:
* FASTERQ.sh: Contains that downloads the raw reads from the SRA
* Preprocess.py: Contains all the functions used in the creation of the final .bam file
* SNV.py: Contains the conversion from the .bam file to the .vcf file
* CNA.py: Contains the conversion from the .bam file to the .csv file
* S2M2.py: Stands for SeCNV to MEDICC2, converts the .csv file produced from SeCNV into the usable .tsv input file for MEDICC2(and other methods that take a .tsv file)
* VCFtoTSV: Takes a .vcf file and converts it to either a binary or ternary tsv file, with the separator being either a tab or a space, with the option to include locations in the file.
## Supplemental Scripts
Some scripts also included but are not used in the pipeline are the following:
* Plot.py: Takes a .csv file and created a heatmap of the copy numbers.
* Run2Lib.py: takes an either an input file or a directory and a mapping file of a:b and replaces every instance of a in the file/directory with b.
* ValidateAll.py: Validates all .bam files in a given directory.
* SCOPE2CSV.py: converts the output of either the single or paired SCOPE Rscript into a .csv file in the same style as one generated by SeCNV.
* RandCells.py: Program that reduces the size of a vcf file using random sampling of SNV cells.
* RandSample.py: Program that reduces the size of a vcf file using random sampling of SNV loci.