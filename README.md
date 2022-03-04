# DeNovoGenomeAssembly
A pipeline for De Novo Genome Assembly

This repository includes this README, which gives a detailed outline of the steps needed perform de novo genome assembly, and the script to a generalized genome assembly pipeline. This pipeline takes in short, paired-end reads from bacteria. A quality assessment is done before and after the assembly to ensure proper contig generations. The resulting contigs can then be used for gene prediction.

**Authors**: Paramita Chatterje, Akul Chopra, Katherine Duchesneau, Xu Qiu, Saideep Narendrula, Haojun Song, Huy Tran, Joseph Tsenum
## Installation with Conda
This assembly uses conda to install and manage the software tools listed below. It is recommended to use this management system to easily ensure the required dependencies for the software tools are properly installed. For more information, visit https://docs.conda.io/en/latest/. 
(Refer to the "References" section for assembly tools documentation and alternative installation methods.)
```bash
conda install -c bioconda fastp multiqc idba megahit spades
```
**Genome assembly tools used:**
- FastP v0.23.1
- MultiQC v1.6
- IDBA v1.1.3
- MEGAHIT v1.2.9
- SPAdes v3.15.4
- QUAST v5.0.2

## 1. Pre-quality Assessment with FastP
Before assembly, we quality checked and filtered samples using FastP. FastP provides an inclusive quality control, it filters out low quality reads, trims reads, cuts adapters, and corrects mismatched base pairs. Additionally, it support multithreading. FastP accepts fastq files, and produces trimmed and filtered fastq files. In addition, it provides a quality report of the reads before and after filtering results, including base content and read quality graphs. 

FastQC requires input fastq files and an output directory. It also allows for extensive customization of the trimming process, a subset of which includes:
- ```-i```: read1 input file name
- ```-I```: read2 input file name
- ```-o```: read1 output file name
- ```-O```: read2 output file name
- ```-M```: mean quality requirement for cutting
- ```-5```: cut bases from front based on mean quality score
- ```-3```: cut bases from tail based on mean quality score
- ```-t```: number of threads
- ```-q```: Phred quality score

```bash
fastp -i <read1.fq> -I <read2.fq> -o <output_name1> -O <output_name2> -M <num> -5 -3 -t <num> -q <num>  
```
This will create multiple filtered fastq files and an html file containing the quality scores.

## 2. *de Novo* Assembly Using IDBA-ud, MEGAHIT, SPAdes
IDBA-UD, MEGAHIT, and SPAdes are *de Novo* assemblers that implement a multiple k-mer strategy to create an iterative de Bruijn graph. Specifically, IDBA-UD and SPAdes creates a basic de Bruijn graph, while MEGAHIT creates a succinct de Bruijn graph (SdBG). The multiple k-mer strategy involves using small k-mer sizes to fill in gaps, while large k-mer sizes resolves any repeats in the sequence. 

#### 2.1. IDBA-UD
Before assembly, IDBA-UD requires paired-end reads to be merged and stored into a single FASTA file. Use fq2fa to merge your reads into a single file.
```bash
fq2fa --merge <read1.fq> <read2.fq> <read.fa>  
```
The resulting FASTA file will then be used for the IDBA-UD input.
- ```-r```: merged FASTA file
- ```-o```: output file name
- ```--mink```: minimum k-mer value
- ```--maxk```: maximum k-mer value
- ```--step```: increment number for increasing k-mer iterations
- ```--min_contig```: minimum contig value
```bash
idba-ud -r <read.fa> -o <output_name> --mink <num> --maxk <num> --step <num> --min_contig <num>
```
IDBA-UD outputs a final contig file, a scaffolds file, and the contig file for every de Bruijn graph it made during the iteration.

#### 2.2. MEGAHIT
MEGAHIT is described as an ultra-fast and memory-efficient NGS assembler. While it is optimized for metagenomes, it also works on generic single genome assembly or single-cell assembly.
- ```-1```: read1 input file name
- ```-2```: read2 input file name
- ```-o```: output file name
```bash
megahit -1 <read1.fq> -2 <read2.fq> -o <output_name> 
```
The output will contain final.contigs.fa, alongside other checkpoints and intermediate data of the SdBG.
#### 2.3. SPAdes
SPAdes was initially designed for small genomes. It was tested on bacterial (both single-cell MDA and standard isolates), fungal, and other small genomes. SPAdes is not intended for larger genomes (e.g., mammalian size genomes).
- ```-k```: list of k-mer lengths
- ```-1```: read1 input file name
- ```-2```: read2 input file name
- ```-o```: output directory location
- ```--isolate```: isolate mode for high-coverage isolates and multi-cell Illumina data (optional)
- ```--careful```: careful mode that reduces number of mismatches and short indels(optional)
```bash
spades.py -k <num,num,...> --<mode> -1 <read1.fq> -2 <read2.fq> -o <output_directory>
```
SPAdes will store all output files in the <output_directory>, of which it will create scaffolds, contigs, assembly graphs, and paths in the assembly graph.
## 3. Post-quality Assessment with QUAST
QUAST is an evaluation tool for assemblies that presents several metrics with or without a reference genome.
- ```-o```: output directory name
- ```-r```: reference genome file (optional)
- ```--min-contig```: minimum threshold for contig length
- ```--features```: genomic feature positions in the reference genome
```bash
quast.py -o <output_directory> -r <reference.fna> -min-contig <num> --features <annotation.gff>
```
This will output report files with tables of summary statistics for each sample isolate.
## References
- [FastP: an ultra-fast all-in-one FASTQ preprocessor](https://github.com/OpenGene/fastp)
- [MultiQC: a summary statistics reporting tool](https://multiqc.info/docs/)
- [IDBA-UD: a basic iterative de Bruijn graph assembler](https://denbi-metagenomics-workshop.readthedocs.io/en/latest/assembly/idba_ud.html)
- [MEGAHIT: an ultra-fast and memory-efficient (meta-)genome assembler](https://github.com/voutcn/megahit)
- [SPAdes: St. Petersburg genome assembler](http://cab.spbu.ru/software/spades/)
- [QUAST: Quality Assessment Tool](http://quast.sourceforge.net/docs/manual.html)
