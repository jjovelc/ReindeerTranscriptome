# Quantify single cell libraries with Salmon-Alevin

## Data availability
All files described here can be found at:
https://uofc-my.sharepoint.com/:f:/r/personal/juan_jovel_ucalgary_ca/Documents/jj/UofC/data/jeffBiernaskie/reindeer_transc_assembly?csf=1&web=1&e=bSfV1Y


Due to the lack of a reindeer reference transcriptome (or a genome with a well-annotated GTF file for cell ranger), the cow (Bos taurus) genome and GTF file have been used before with Cell Ranger for quantification of 10X libraries.

Our assembled reindeer transcriptome was aligned against the cow transcriptome and those putative novel reindeer transcripts that did not align against that were deemed novel transcripts and were contatenated to the cow transcriptome to constitute the hybrid `bt-rt_transcriptome.fa`. 

When using a custom transcriptome as reference to quantify single cell libraries, Salmon is a good choice because it is flexible to create indices and fast to run quantification with its module 'alevin'. Here we describe the procedure to use the hybrid reindeer transcriptome to quantify single cell libraries.

## Index transcriptome

For quantification, `salmon alevin` needs two pieces of information.

1. An index of the reference transcriptome
2. A transcript-to-gene list (t2g.txt)

### Index reference transcriptome

At the moment of writing this repo, the last version of Salmon 1.10.1 fixed a bug in quantification in previous versions. However, such version is not yet available in Anaconda. For that reason an aptainer (singularity) container was created for such a version, since installing from source code in ARC is often troublesome.

The container can be run interactively with the following command:

```bash
    singularity run --cleanenv salmon_singularity.sif
```
However, it is more practical to submit to an ARC compute node using the slurm file: `format_index_salmon.slurm`. If you want to format a different transcriptome change inside the slurm file the input file `bt-rt_transcriptome.fa` accordinginly. Make sure your reference transcriptome file has the extension 'fa' or 'fasta' or 'fna'.

The index directory will be `bt-rt_transcriptome_idx` and contains all the files needed for `salmon alevin` for quantification. 

### Create transcript-to-gene list (t2g.txt)

For that, simply run the following command:

```bash
    perl create_t2g.pl bt-rt_transcriptome.fa > t2g.txt
```

## Quantify libraries with salmon alevin

To quantify a series of libraries with barcodes-UMIs in end1 and reads in end2, you can run the following command:

```bash
    sbatch quantify_with_salmon.slurm
```

Results from `salmon alevin` will be stored in directory `alevin_output`. Now move this directory to your local computer.

To inspect QC metricts on your alevin results, you can use, in your local computer, the script alevinQC.R. Run it in R studio, since running from the prompt may result in error when generating the HTML report. You only have to modify the working directory to include the absolute path up to your alevin_output directory. You then can modify the value of the argument `sampleId` to any string you want to use to identify this run. 

Finally, to import data into R, for example with use with Seurat, do something like this:

```R
library(Seurat)
library(Matrix)
library(data.table)
library(DropletUtils)

setwd('/your/working/directory/alevin_output/')

# Path to the Alevin output directory
alevin_dir <- "alevin"

# File paths
files <- file.path(alevin_dir, "quants_mat.gz")
names(files) <- "yourSampleName"

# Import with tximport
txi <- tximport(files, type = "alevin")
data <- txi$counts

# Create the Seurat object
pbmc <- CreateSeuratObject(counts = data, project = "YourProjectName")
```