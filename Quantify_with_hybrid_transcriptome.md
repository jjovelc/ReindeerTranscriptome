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

