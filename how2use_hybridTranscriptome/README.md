# Using the hybrid transcriptome

## Introduction
  Initially, a series of transcriptomes assemblies were conducted as described [here](../README.md). In a nutshell, several transcriptome assemblies were conducted using three datasets. One of those assemblies was conducting using data from bulk RNAseq libraries from skin and antlers of reindeer, as described in [Sinna et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36493752/). Because the main goal of assembling this hybrid (*Rangifer tarandus* & *Bos taurus*) was to improve on the *B. taurus* transcriptome being used as reference for RNAseq studies in *R. tarandus*, we gave 'preference' to those transcripts assembled from *R. tarandus* data from skin and antlers (Cell paper referenced above). Below as schematics of the assembly process.

  <p align="center">
  <img src="flowDiagram.png" alt="Hybrid assembly pipeline" title="Assembly pipeline" width="50%">
</p>

Briefly, we conducted a combined assembly that included data from the three sources described [here](../README.md). In parallel, we conducted another assembly with only   
Jeff's data (skin and antler libraries). We then aligned both, combined and Jeff, assemblies against the *B. taurus* transcriptome. Only those transcripts that did not align to such reference were used for further processing. The hybrid transcriptome was assembled by the following sets of sequences:

- 1. The entire Ensembl *B. taurus* transcriptome (Bos_taurus.ARS-UCD1.2.cdna.all)
- 2. Jeff's transcripts that did not aligned against (i)
- 3. Transcripts in the combined assembly that did not align against (i) and (ii)

## How to use hybrid transcriptome

Just as a reminder, singularity container for salmon 1.10.1 was created because a quantification bug was reported in the version that is distributed through conda. Once this problem is fixed, the conda version should be OK to use.

The container containing the salmon software can be found here:
https://uofc-my.sharepoint.com/:u:/r/personal/juan_jovel_ucalgary_ca/Documents/jj/UofC/data/jeffBiernaskie/reindeer_transc_assembly/salmon_singularity.sif?csf=1&web=1&e=5X0QjN

The hybrid transcriptome can be found here:
https://uofc-my.sharepoint.com/:u:/r/personal/juan_jovel_ucalgary_ca/Documents/jj/UofC/data/jeffBiernaskie/reindeer_transc_assembly/bt-rt_transcriptome.fa?csf=1&web=1&e=8c5bW9

Make sure that the container and the transcriptome files are in the same folder.

1. Generate a transcript-to-gene list required for salmon alevin, with the following command: 
   ```bash
      perl create_transc2gene.pl bt-rt_transcriptome.fa > txp2gene.tsv
   ```

   Relevant script is here: [create_transc2gene.pl](../create_transc2gene.pl)

2. Format the transcriptome using script [format_index_salmon.slurm](../format_index_salmon.slurm), simply by running command:

```bash
  sbatch format_index_salmon.slurm
``` 

3. Quantify scRNAseq libraries using script [quantify_with_salmon.slurm](../quantify_with_salmon.slurm), simply by running command:

```bash
  sbatch quantify_with_salmon.slurm
```

This script assumes that the libraries are in PWD and their names end in R1/2.fq.gz; if that is not your case, modify the following loop in the script:

```bash
  for BARCODES in ${DIR}/*_R1.fq.gz; do
    BARCODES_FILES+=($BARCODES)
    READS="${BARCODES%_R1.fq.gz}_R2.fq.gz"
    READS_FILES+=($READS)
  done
```

This will generate, an output directory per sample. 

In a pilot test, when comparing the quantification of 6 scRNAseq libraries (provided by Ross), with either the *B. taurus* or the hybrid transcriptome, the following amount of genes were detected in each case:

|  Transcriptome  | Sample      | # genes detected |
|-----------------|-------------| -----------------|
| *B.taurus*      | R_BM_1      |    13,980        |  
|                 | R_BM_2      |    14,244        |
|                 | R_BM_A      |    14,238        |
|                 | R_BM_S      |    14,236        |
|                 | Skin_R_N    |    16,469        |
|                 | Skin_R_PM   |    15,342        |
|                 |             |                  |
| Hybrid          | R_BM_1      |    15,884        |
|                 | R_BM_2      |    16,209        |
|                 | R_BM_A      |    16,205        |
|                 | R_BM_S      |    16,204        |
|                 | Skin_R_N    |    18,880        |
|                 | Skin_R_PM   |    17,468        |


Those directories can be moved to a local computer and further processed in R. For an example of how to import and further process quantified libraries in R, see script [BM_skin_postprocessing.R](../miscelaneous_scripts/BM_skin_postprocessing.R).

If you have questions or suggestions, please contact me (juan.jovel@ucalgary.ca).

