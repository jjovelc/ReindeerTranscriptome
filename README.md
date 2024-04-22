# ReindeerTranscriptome

### Principal Investigator: Jeff Biernaskie
### Institution: Faculty of Veterinary Medicine, University of Calgary

### Important notes
All data mentioned in this document can be found at:
https://uofc-my.sharepoint.com/:f:/r/personal/juan_jovel_ucalgary_ca/Documents/jj/UofC/data/jeffBiernaskie/reindeer_transc_assembly?csf=1&web=1&e=SdzCaP

All scripts mentioned here and not included in full are stored in this repo.

## Description of Project

In adult mammals, deep skin injuries heal by rapid repopulation with reactive fibroblasts that deposit extracellular matrix and form fibrotic scar tissue. Discovery of a tight-skinned, free-living mammal that exhibits both skin regeneration and fibrotic repair would enable comparative insights that better contextualizes human healing. In reindeer (Rangifer tarandus), antlers are regenerated annually in both sexes and grow at an explosive rate exceeding 1 cm in length each day. Growing antlers are covered by specialized skin called velvet. Therefore, antler velvet itself may harbor innate regenerative capacity and could represent a unique model to study the molecular events enabling adult skin regeneration.

Despite considerable endeavors, a well annotated genome and transcriptome are still missing for the reindeer. In collaboration with Drs. Claude Robert and Julien Prunier (University of Laval) and Dr. Juha Kantanen (Natural Resources Institute of Findland), we initiated an effort to assembly the reindeer transcriptome. Here we describe the bioinformatics procedures undertaken.

## Bioinformatics

### Description of datasets

The main directory of the project in ARC is:
'/work/vetmed_data/jj/projects/jeffBiernaski/reindeer/transcriptome_assembly/'

1. Jeff Biernaskie data: 35 single-end (150 bp) RNAseq libraries from antlers and back tissue, described in (https://pubmed.ncbi.nlm.nih.gov/36493752/).
2. Julien Prunier data: 1 paired-end (150 bp) RNAseq libraries from organ tissues.
3. Juha Kantanen data: 46 paired-end (75 bp) RNAseq libraries from 3 adipose tissues of a male.

### Preliminary assemblies

Three well-reputated assembly software were tested to assemble each of the transcriptomes:

#### Trinity
Model: The Bruijn graph-based
Desirable features: Downstream processing and annotation tools
Developers: Team of Aviv Regev and Brian Hass (MIT)
Reference: Nat. Protoc. 8, 1494-1512. (2013)

#### RNAspades
Model: The Bruijn graph and graph-based 
Desirable features: Very successful assemblers for genome and metagenome assembly
Developers: Team of Andrey D Prjibelski (Center for Algorithmic Biotechnology)
Reference: Gigascience 8:9, giz100 (2019)

#### Transabyss
Model: The Bruijn graph and graph-based
Desirable features: Laval group had good experience with it. 
Developers: Team of Inanc Birol (BC Cancer Research Institute)
(Inanc Birol) Nat. Methods 7, 909-912. (2010).

Assemblies were conducted separately, for each data set, and are stored in the following subdirectories:

jeffData <br>
juhaData <br>
julienData <br>

##### Quality trimming

```bash
# For Jeff data:
fastq-mcf ~/useful_files/adapters.fa -o rangifer_tarandus_RNAseq_jeff_trim30.fq rangifer_tarandus_RNAseq_jeff.fastq -k 0 -l 50  -w 3 -q 30

# For Juha data:
fastq-mcf ~/useful_files/adapters.fa -o rangifer_tarandus_RNAseq_juha_trim35_R1.fq -o rangifer_tarandus_RNAseq_juha_trim35_R2.fq rangifer_tarandus_RNAseq_juha_R1.fq  rangifer_tarandus_RNAseq_juha_R2.fq -k 0 -l 50  -w 3 -q 35

# For Julien data:
fastq-mcf ~/useful_files/adapters.fa -o rangifer_tarandus_RNAseq_julien_trim30_R1.fq -o rangifer_tarandus_RNAseq_julien_trim30_R2.fq rangifer_tarandus_RNAseq_julien_R1.fq  rangifer_tarandus_RNAseq_julien_R2.fq -k 0 -l 50  -w 3 -q 30
```

##### Assembly

```bash
# RNAspades

rnaspades.py -t 48  -s rangifer_tarandus_RNAseq_jeff_trim30.fq -o Rt_jeff_rnaspades_assembly

declare -a infiles=("rangifer_tarandus_RNAseq_juha_trim35_R1.fq" "rangifer_tarandus_RNAseq_julien_trim30_R1.fq")

for file in "${infiles[@]}"
do
    rnaspades.py -t 48 -1 "$file" -2 "${file/_R1/_R2}" -o "$(basename "$file" _R1.fq)"_rnaspades_assembly
done


# TransAbyss

transabyss --threads 48 --se rangifer_tarandus_RNAseq_jeff_trim30.fq

transabyss --threads 48  --pe rangifer_tarandus_RNAseq_juha_trim35_R1.fq rangifer_tarandus_RNAseq_juha_trim35_R2.fq

transabyss --threads 48  --pe rangifer_tarandus_RNAseq_julien_trim30_R1.fq rangifer_tarandus_RNAseq_julien_trim30_R2.fq

# Trinity

Trinity --seqType fq  --single rangifer_tarandus_RNAseq_jeff_trim30.fq --CPU 48 --max_memory 256G --trimmomatic --output Rt_jeff_trinity_assembly

declare -a infiles=("rangifer_tarandus_RNAseq_juha_trim35_R1.fq" "rangifer_tarandus_RNAseq_julien_trim30_R1.fq")

for file in "${infiles[@]}"
do
    rnaspades.py -t 48 -1 "$file" -2 "${file/_R1/_R2}" -o "$(basename "$file" _R1.fq)"_rnaspades_assembly

    Trinity --seqType fq --left "$file" --right "${file/_R1/_R2}" --CPU 48 --max_memory 256G --trimmomatic --output "$(basename "$file" _R1.fq)"_trinity_assembly
done

```

#### Evaluation of assembly

```bash
# QUAST

declare -a names=("jeff" "juha" "julien")

for name in "${names[@]}"
do
    quast.py transabyss_2.0.1_assembly/transabyss-final.fa -o "Rt_${name}_transabyss_quast" --fast
done

# BUSCO
declare -a NAMES=("jeff" "juha" "julien")
# A directory with symbolic links was created
DIR="/work/vetmed_data/jj/projects/jeffBiernaski/reindeer/transcriptome_assembly/simb_links/"
LINNEAGE="${DIR}/cetartiodactyla_odb10"

for NAME in "${NAMES[@]}"
do
    busco -i "${DIR}/${NAME}_transabyss-final.fa" -l "$LINNEAGE" -o "Rt_${NAME}_transabyss_busco" -m transcriptome -c 48
done

```

Here a summary of assembly results:

![Assembly results](summaryAssemblies.png "Assembly results")


From the statistics presented above, it is clear that RNAspades and Trinity performed better on all datasets. We had the additional complication that each set of libraries are different: paired-end or single-end, different length and different depth.

We could not make any of the aligners work with a simultaneously using single-end and paired-end reads. Since Jeff's data is single-end, we decided to go with single-end data for all datasets and try to downsize each dataset to the size of the smallest one.

### Single-end assembly

Data in ARC:

'/work/vetmed_data/jj/projects/jeffBiernaski/reindeer/transcriptome_assembly/all_sequences_together'

In order to improve quality of assembly, we trimmed data using a Q threshold of 40.

```bash
# Jeff data:
fastq-mcf /home/juan.jovel/useful_files/adapters.fa jeff_raw_data_all_trim35.fq -o jeff_raw_data_all_trim40.fq  -l 75 -q 40

# Juha data:
fastq-mcf /home/juan.jovel/useful_files/adapters.fa juha_raw_data_all_1_trim35.fq juha_raw_data_all_2_trim35.fq -o juha_raw_data_all_1_trim40.fq -o juha_raw_data_all_2_trim40.fq  -l 50 -q 40

# Julien data:
fastq-mcf /home/juan.jovel/useful_files/adapters.fa julien_raw_data_all_1_trim35.fq julien_raw_data_all_2_trim35.fq -o julien_raw_data_all_1_trim40.fq -o julien_raw_data_all_2_trim40.fq  -l 100 -q 40
```

For paired-end data, end1 and end2 were merged as follows:

```bash
declare -a INFILES=("julien_raw_data_all_1_trim40.fq" "juha_raw_data_all_1_trim40.fq")

for FILE in "${INFILES[@]}"
do
    bbmerge-auto.sh in1="$FILE" in2="${FILE/_1/_2}" out="${FILE/_1*/}_merged_trim40.fq" adapters=/home/juan.jovel/useful_files/adapters.fa
done
```

After rezising files we ended up withthe following files:

52G  jeff_raw_data_all_trim40_150M.fq
58G  juha_raw_data_all_merged_trim40_250M.fq
54G  julien_raw_data_all_merged_trim40.fq

Initially, assembly was attempted with RNAspades, but it finished without any error but also without reporting any assembly. We then moved to attempt the assembly with Trinity as follows:

```bash
Trinity --seqType fq --single jeff_raw_data_all_trim40_150M.fq,juha_raw_data_all_merged_trim40_250M.fq,julien_raw_data_all_merged_trim40.fq --CPU 48 --max_memory 250G --trimmomatic --output 60Gb_Rt_all_trinity_single_assembly
```

Trinity was able to assemble 898,563 longer than 200 bp (file `Trinity.fasta`). We then applied TransDecoder to initially scan for long open reading frames. `TransDecoder.LongOrfs` was run with the following command:

```bash
    TransDecoder.LongOrfs -t transcripts.fasta
```

This tool works by scanning transcript sequences for ORFs that are long enough to be potential protein-coding regions. It typically looks for standard start and stop codons and uses the length of the ORF as one of the criteria for predicting coding regions.

In a second step, we applied the sister program `TransDecoder.Predict`. TransDecoder.Predict is designed to analyze the ORFs identified by TransDecoder.LongOrfs (or other methods) and determine which of these are most likely to represent actual protein-coding regions.  TransDecoder.Predict searches for the presence of known protein domains.

```bash
   TransDecoder.Predict -t transcripts.fasta  
```

TransDecoder.Predict implicitly relies on the intermediate files generated by TransDecoder.LongOrfs. It uses these files to analyze the ORFs and predict which ones are likely to be true protein-coding regions. By default, it looks at directory `Trinity.fasta.transdecoder_dir` in the same directory where TransDecoder.LongOrfs was run and uses the intermediate files generated by this latter one to predict putative coding ORFs.

This resulted in 152,529 predicted proteins, which are stored in files `Trinity.fasta.transdecoder.cds` and `Trinity.fasta.transdecoder.pep` containing, respectively, DNA and protein predicted sequences.

From those, only complete ORFs were extracted (5' and 3' moeities were excluded). This produced 80,034 putatively novel transcripts, and they are stored in file `Trinity.fasta.transdecoder_complete.cds` and `Trinity.fasta.transdecoder_complete.pep`.

Complete ORfs were extracted with the following commands:

```bash
    perl extract_complete_ORFs.pl Trinity.fasta.transdecoder.cds > Trinity.fasta.transdecoder_complete.cds

    perl extract_complete_ORFs.pl Trinity.fasta.transdecoder.pep > Trinity.fasta.transdecoder_complete.pep

```

More generally, if any subset of a FASTA file needed to be extracted, the subset of IDs was extracted and used for pulling out from a FASTA file as follows:

```bash
    # For nucleotide sequences
    perl extract_nt_records.pl <fasta-file> <ids-file> > <output.fasta>

    # For protein sequences
    perl extract_protein_records.pl <fasta-file> <ids-file> > <output.fasta>
```

Such putative transcripts were subjected to annotation with the Trinotate pipeline (https://github.com/Trinotate/Trinotate/wiki). 

Initially, an sqlite database is created and populated with databases for comparison of nucleotide and protein sequences. 

```bash
Trinotate --db reindeer_transc.db --create --trinotate_data_dir trinotate_db
```

Then, sqlite database was initialized by importing nucleotide and protein sequences.

```bash
Trinotate --db reindeer_transc.db \
        --init \
        --gene_trans_map "Trinity.fasta.transdecoder_complete.cds.gene-map.txt" \
        --transcript_fasta "Trinity.fasta.transdecoder_complete.cds" \
        --transdecoder_pep "Trinity.fasta.transdecoder_complete.pep"
```

 The gene-map file was created with the following command:

```bash
# generation of a gene-map file
awk '{ if ($0 ~ /^>/) { id=substr($0,2); print id "\t" id } }'  Trinity.fasta.transdecoder_complete.cds > Trinity.fasta.transdecoder_complete.cds.gene-map.txt

```

Trinotate was run:

```bash
DB_PATH="reindeer_transc.db"

Trinotate --db ${DB_PATH} \
--run ALL \
--CPU 12 \
--transcript_fasta "Trinity.fasta.transdecoder_complete.cds" \
--transdecoder_pep "Trinity.fasta.transdecoder_complete.pe"

Trinotate --db ${DB_PATH} --LOAD_swissprot_blastp uniprot_sprot.ncbi.blastp.outfmt6
Trinotate --db ${DB_PATH} --LOAD_pfam TrinotatePFAM.out
Trinotate --db ${DB_PATH} --LOAD_signalp6 sigP6outdir/prediction_results.txt
Trinotate --db ${DB_PATH} --LOAD_EggnogMapper eggnog_mapper.emapper.annotations
Trinotate --db ${DB_PATH} --LOAD_tmhmmv2 tmhmm.v2.out

Trinotate --db ${DB_PATH} --LOAD_swissprot_blastx uniprot_sprot.ncbi.blastx.outfmt6
Trinotate --db ${DB_PATH} --LOAD_infernal infernal.out
```
And results were extracted into a tabular report:

```bash
    Trinotate --db  reindeer_transc.db --report > reindeer_transc_Trinotate_report.tsv
```

Reads for which Trinotate found annotations can be extracted with the following command line:

```bash
    cut -f 1 Trinity.fasta.transdecoder_complete.cds | sed 's/>//' | grep - reindeer_transc_Trinotate_report.tsv | cut -f 1,3 | awk '$2 != "."' > Trinity.fasta.transdecoder_complete_trinotate_annotations.txt
```

Out of 80,034 reads subjected to annotation in file `Trinity.fasta.transdecoder_complete.cds`, 43,694 got annotations and those annotations were saved in file `Trinity.fasta.transdecoder_complete_annotations.txt`.  

Unique protein IDs were extracted from the Trinotate report with the following command:

```bash
cut -f 3 Trinity.fasta.transdecoder_complete_annotations.txt | cut -d '^' -f 1 | sed '/^\.$/d' | sort | uniq > unique_reindeer_protein_ids_list.txt
```

Annotation for such proteins were retrieved from the Uniprot database:

```bash
python retrieve_protein_annotations_uniprot.py unique_reindeer_protein_ids_list.txt > reindeer_assembly_annotations_uniprot.tsv 
```

The following files were used for annotation the ID of files:

reindeer_assembly_annotations_uniprot.tsv
Trinity.fasta.transdecoder_complete_annotations.txt
node_and_prot_ids.txt
Trinity.fasta.transdecoder_complete.cds
Trinity.fasta.transdecoder_complete.pep

using script `annotate_fasta_id.pl`. The results were stored in the following files:

trinotate_report_table.tsv
trinotate_report_table.xlsx 
(a table with annotations for each transcript that could be annotated)

Trinity.fasta.transdecoder_complete_with_annotation.fasta
(All nuleotide reads that could be annotated)

Trinity.fasta.transdecoder_complete_with_annotation.faa
(All protein reads that could be annotated)

Annotated sequences were aligned agains the `Bos_taurus.ARS-UCD1.2.cdna.all.fa` transcriptome and transcripts that did not aligned that reference were labelled as novel transcripts. Those are stored in file: `Trinity.fasta.transdecoder_complete_with_annotation_novel.fasta`.

However, it is important to notice that even when those transcripts aligned to the cow transcriptome with some similarity and significant e-value it does noe mean that the could not align with higher similarity to the reindeer actual transcriptome.
 








