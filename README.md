---

#Ajay Kumar Mahato#




---

---
#Genome Assembly Using Soapdenovo#
---

**Genome Assembly:**

In bioinformatics, sequence assembly refers to aligning and merging fragments from a longer DNA sequence in order to reconstruct the original sequence. When the genome of a species is to be sequenced, the chromosomes from many cells are broken at random positions into small fragments, which are sequenced, and reassembled into long sequences (contigs). Contigs may be assembled into longer sequences called scaffolds and sometimes, if the depth of sequencing is high enough, there may be enough information to assemble most of the scaffolds into chromosomes. The resulting collection of sequences after assembly is called a genome assembly
---

**Contig:**

A contig is a contiguous stretch of DNA sequence without gaps that has been assembled solely based on direct sequencing information. Short sequences (reads) from a fragmented genome are compared against one another, and overlapping reads are merged to produce one long sequence. This merging process is iterative: overlapping reads are added to the merged sequence whenever possible and so the merged sequence becomes even longer.When no further reads overlap the long merged sequence, then this sequence called a contig has reached its maximum length.
---


**Scaffold:**


Supercontigs or scaffolds are sets of ordered, oriented contigs which are joined together on the basis of PE/MP reads information and the gaps are represented by variable number of Ns.

---

#List of software widely used for assembly#

```bash
    A) Genome = MIRA, VELVET, SOAPdenovo, ABySS, SPADES, FALCON, CANU etc.

    B) Transcriptome= Trinity, Trans-ABySS, VELVET-Oases, etc.
```
---

**Practicle Session**
---

##Genome Assembly using SOAPdenovo##

SOAPdenovo aims for large plant and animal genomes, although it also works well on bacteria and fungi genomes.

	1: It runs on 64-bit Linux system with a minimum of 5G physical memory. For big genomes like human, about 150 GB memory would be required.
    2: Source code of 63mer and 127mer versions were merged. But two versions of executable program are provided. 
    3: The 63mer version support kmer only ≤63. 
    4: The 127mer version support kmer only ≤127 and double the memory consumption than 63mer version, even being used with kmer ≤63.

---

Step by Step Tutorial
---



Step 1: Install SOAPdenovo using command line
---
```bash
ajay@ajay-550P5C-550P7C:~$ sudo apt-get install soapdenovo
[sudo] password for ajay:        
Reading package lists... Done
Building dependency tree       
Reading state information... Done
The following NEW packages will be installed:
  soapdenovo
0 upgraded, 1 newly installed, 0 to remove and 54 not upgraded.
1 not fully installed or removed.
Need to get 347 kB/14.2 MB of archives.
After this operation, 977 kB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu bionic/universe amd64 soapdenovo amd64 1.05-4 [347 kB]
Fetched 347 kB in 1s (235 kB/s)     
Selecting previously unselected package soapdenovo.
(Reading database ... 318779 files and directories currently installed.)
Preparing to unpack .../soapdenovo_1.05-4_amd64.deb ...
Unpacking soapdenovo (1.05-4) ...
Setting up soapdenovo (1.05-4) ...
Processing triggers for man-db (2.8.3-2ubuntu0.1) ...
ajay@ajay-550P5C-550P7C:~$
```
---


Step 2- prepare configuration file (configuration file stores all the parameters used for running the assembly like; data type, library insert size, read length etc.)
---


```bash

 SOAPdenovo full Configuration file details

#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=200
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=/path/**LIBNAMEA**/fastq1_read_1.fq
q2=/path/**LIBNAMEA**/fastq1_read_2.fq
#another pair of fastq file, read 1 file should always be followed by read 2 file
q1=/path/**LIBNAMEA**/fastq2_read_1.fq
q2=/path/**LIBNAMEA**/fastq2_read_2.fq
#a pair of fasta file, read 1 file should always be followed by read 2 file
f1=/path/**LIBNAMEA**/fasta1_read_1.fa
f2=/path/**LIBNAMEA**/fasta1_read_2.fa
#another pair of fasta file, read 1 file should always be followed by read 2 file
f1=/path/**LIBNAMEA**/fasta2_read_1.fa
f2=/path/**LIBNAMEA**/fasta2_read_2.fa
#fastq file for single reads
q=/path/**LIBNAMEA**/fastq1_read_single.fq
#another fastq file for single reads
q=/path/**LIBNAMEA**/fastq2_read_single.fq
#fasta file for single reads
f=/path/**LIBNAMEA**/fasta1_read_single.fa
#another fasta file for single reads
f=/path/**LIBNAMEA**/fasta2_read_single.fa
#a single fasta file for paired reads
p=/path/**LIBNAMEA**/pairs1_in_one_file.fa
#another single fasta file for paired reads
p=/path/**LIBNAMEA**/pairs2_in_one_file.fa
#bam file for single or paired reads, reads 1 in paired reads file should always be followed by reads 2
#   NOTE: If a read in bam file fails platform/vendor quality checks(the flag field 0x0200 is set), itself and it's paired
read would be ignored.
b=/path/**LIBNAMEA**/reads1_in_file.bam
#another bam file for single or paired reads
b=/path/**LIBNAMEA**/reads2_in_file.bam
[LIB]
avg_ins=2000
reverse_seq=1
asm_flags=2
rank=2
# cutoff of pair number for a reliable connection (at least 5 for large insert size)
pair_num_cutoff=5
#minimum aligned length to contigs for a reliable read location (at least 35 for large insert size)
map_len=35
q1=/path/**LIBNAMEB**/fastq_read_1.fq
q2=/path/**LIBNAMEB**/fastq_read_2.fq
f1=/path/**LIBNAMEA**/fasta_read_1.fa
f2=/path/**LIBNAMEA**/fasta_read_2.fa
p=/path/**LIBNAMEA**/pairs_in_one_file.fa
b=/path/**LIBNAMEA**/reads_in_file.bam
```

__Notable fields include average insert size and read length, which differ depending on the sequencing technology, and q1, q2, and q; the paths to the forward, reverse and singles trimmed reads.__


---
Step 3: Format of configuration file used during this exercise
---

```bash
#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=350
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 250 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
# path to genes
q1=/PATH TO YOUR FOLDER WHICH CONTAINS THE DOWNLOADED ILlumina FORWADED READS/SRR1928200_1.fastq
q2=/PATH TO YOUR FOLDER WHICH CONTAINS THE DOWNLOADED ILlumina reverse READS/SRR1928200_2.fastq

```
---
Running SOAPdenovo for de-novo assembly
---

To run the assembler we will use the SOAPdenovo-63mer command with the all option (to perform kmer graph construction, contig error correction, mapping of reads to contigs, and scaffolding)

```bash

-s = for the path to the config file
-K=  for the size of the kmer
-o =  for the output prefix
1=  for assembly log,
2= for assembly errors.

```
---

Command for executing SOAPdenovo
---
```bash

$ SOAPdenovo-63mer all -s .config -K 31 -R -o graph_Sample_31 1>ass31.log 2>ass31.err

```
---
# repeat genome assembly  for k=35, k=41, etc to get optimal genome assembly.
