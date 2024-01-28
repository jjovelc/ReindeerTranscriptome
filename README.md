# ReindeerTranscriptome

### Principal Investigator: Jeff Biernaskie
### Institution: Faculty of Veterinary Medicine, University of Calgary

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




