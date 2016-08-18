[![Github Issues](http://githubbadges.herokuapp.com/bretonics/Bioinformatics-Software/issues.svg)](https://github.com/bretonics/Bioinformatics-Software/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/Bioinformatics-Software/pulls.svg)](https://github.com/bretonics/Bioinformatics-Software/pulls)
![](https://reposs.herokuapp.com/?path=bretonics/Bioinformatics-Software&color=orange)


#Workflow Tools
| Program        | Description    | Source           |
| :------------- | :------------- | :-------------   |
| [Artemis](http://www.sanger.ac.uk/science/tools/artemis) | A genome browser and annotation tool that allows visualisation of sequence features, next generation data and the results of analyses within the context of the sequence, and also its six-frame translation. | [Download](http://www.sanger.ac.uk/science/tools/artemis)
| [BamTools](https://github.com/pezmaster31/bamtools) | C++ API & command-line toolkit for working with BAM (Binary SAM file) data. Provides a programmer's API and an end-user's toolkit for handling BAM files. | [Clone](https://github.com/pezmaster31/bamtools/wiki/Building-and-installing)
|[BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) | Command line application suite of BLAST tools that utilizes the NCBI C++ Toolkit. | [Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
| [EDirect](http://www.ncbi.nlm.nih.gov/books/NBK179288/) | An advanced method for accessing the NCBI's set of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a UNIX terminal. | N/A
| [E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25497/) | Entrez Programming Utilities (E-utilities) are a set of nine server-side programs that provide a stable interface into the Entrez query and database system at the NCBI. | N/A
| [Samtools](http://www.htslib.org/) | A suite of programs for interacting with high-throughput sequencing data (HTS) from next generation sequencing data. It consists of three separate repositories:<br><br>**Samtools:** Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format.<br><br>**BCFtools:** Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants.<br><br>**HTSlib:** A C library for reading/writing high-throughput sequencing data. | [Download](http://www.htslib.org/download/)
| [Seqtk](https://github.com/lh3/seqtk) | Fast and lightweight tool for processing sequences in the FASTA or FASTQ format. | [Clone](https://github.com/lh3/seqtk.git)
| [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK242621/) | The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives. | [Download](http://ncbi.github.io/sra-tools/)
| [VCFtools](https://vcftools.github.io/index.html) | Package designed for working with complex genetic variation data in the form of VCF files.| [Download](https://vcftools.github.io/downloads.html)


#Analysis
| Program        | Description    | Purpose          | Source           |
| :------------- | :------------- | :-------------   | :-------------   |
| [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) | *ab initio*, trainable gene prediction in eukaryotic genomic sequences.  | Gene Prediction | [Download](http://bioinf.uni-greifswald.de/augustus/downloads/)
| [BUSCO](http://busco.ezlab.org/) | Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs <br><br>BUSCO provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness based on evolutionarily informed expectations of gene content from near-universal single-copy orthologs selected from [OrthoDB](http://orthodb.org/). | Assembly Quality Assesment | [Download](http://busco.ezlab.org/)
| [Circlator](http://sanger-pathogens.github.io/circlator/) |  Predict and automate assembly circularization and produce accurate linear representations of circular sequences. | Circularize Genome | [Download](https://github.com/sanger-pathogens/circlator/releases/latest)
| [Clustal](http://www.clustal.org/) | Fast and scalable multiple sequence alignment (can align hundreds of thousands of sequences in hours)| MSA | [Download](http://www.clustal.org/omega/#Download)
| [Galaxy](https://galaxyproject.org/) | Web portal for accessible, reproducible, and transparent computational research. | Analysis package | [Download](https://wiki.galaxyproject.org/Admin/GetGalaxy)
| [HMMER](http://hmmer.org/) | Search sequence databases for sequence homologs, and for making sequence alignments, analyzed by using profile hidden Markov models | Detect Homologs | [Download](http://hmmer.org/download.html)
| [Mauve](http://darlinglab.org/mauve/mauve.html) | A system for constructing multiple genome alignments in the presence of large-scale evolutionary events such as rearrangement and inversion. | Genome Aligner| [Download](http://darlinglab.org/mauve/download.html)
| [Mothur](http://www.mothur.org/) |  Expandable software to fill the bioinformatics needs of the microbial ecology community. | Microbial Ecology Pipeline | [Download](http://www.mothur.org/wiki/Download_mothur)
| [MUMmer Package](http://mummer.sourceforge.net/)| Ultra-fast alignment of large-scale DNA and protein sequences. A system for rapidly aligning entire genomes, whether in complete or draft form. <br><br>**MUMmer** is a suffix tree algorithm designed to find maximal exact matches of some minimum length between two input sequences.<br><br>**NUCmer** is a standard DNA sequence alignment. It is a robust pipeline that allows for multiple reference and multiple query sequences to be aligned in a many vs. many fashion.<br><br>**PROmer** is like NUCmer with one exception - all matching and alignment routines are performed on the six frame amino acid translation of the DNA input sequence. | Genome Aligner | [Download](https://sourceforge.net/projects/mummer/files/latest/download?source=files)
| [MUSCLE](http://www.drive5.com/muscle/) | MUSCLE can align hundreds of sequences in seconds. | MSA | [Download](http://www.drive5.com/muscle/downloads.htm)
| [Picard](http://broadinstitute.github.io/picard/) | Set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF. | HTS Toolkit | [Download](https://github.com/broadinstitute/picard/releases/latest)
| [QIIME](http://qiime.org/) | Bioinformatics pipeline for performing microbiome analysis from raw DNA sequencing data. QIIME is designed to take users from raw sequencing data generated on the Illumina or other platforms through publication quality graphics and statistics. | Microbial Ecology Pipeline | [Install](http://qiime.org/install/index.html)
| [QUAST](http://bioinf.spbau.ru/quast) | Evaluates genome assemblies. | Evaluate Genome Assemblies | [Download](https://sourceforge.net/projects/quast/files/latest/download?source=files)
| [T-Coffee](http://www.tcoffee.org/) | A multiple sequence alignment package that can align sequences (Protein, DNA, and RNA) or combine the output of your favorite alignment methods (Clustal, Mafft, Probcons, Muscle...) into one unique alignment (**M-Coffee**). It is also able to combine sequence information with protein structural information (**3D-Coffee/Expresso**), profile information (**PSI-Coffee**) or RNA secondary structures. | MSA | [Download](http://www.tcoffee.org/Projects/tcoffee/#DOWNLOAD)


##PacBio Sequencing
| Program        | Description    | Purpose          | Source           |
| :------------- | :------------- | :-------------   | :-------------   |
| [BLASR](https://github.com/PacificBiosciences/blasr) | PacBio® long read aligner | Sequence Aligner | [Download](https://github.com/PacificBiosciences/blasr/releases/latest)
| [Canu](https://github.com/marbl/canu) | Fork of the **Celera Assembler** designed for high-noise single-molecule sequencing (such as the PacBio RSII or Oxford Nanopore MinION). | Genome Assembly | [Download](https://github.com/marbl/canu/releases/latest)
| [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page)| Celera Assembler is a de novo whole-genome shotgun (WGS) DNA sequence assembler, and can use any combination of platform reads. | Genome Assembly | [Download](https://sourceforge.net/projects/wgs-assembler/files/latest/download?source=files)
| [Cerulean](https://sourceforge.net/projects/ceruleanassembler/) | Cerulean extends contigs assembled using short read datasets like Illumina paired-end reads using long reads like PacBio RS long reads. | Hybrid Assembly| [Download](https://sourceforge.net/projects/ceruleanassembler/files/latest/download)
| [PBSuite](https://sourceforge.net/projects/pb-jelly/) | **PBJelly** is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes. <br><br>**PBHoney** is an implementation of two variant-identification approaches designed to exploit the high mappability of long reads (i.e., greater than 10,000 bp). PBHoney considers both intra-read discordance and soft-clipped tails of long reads to identify structural variants. | Reference Mapping<br><br>Variant Calling| [Download](https://sourceforge.net/projects/pb-jelly/files/latest/download)
| [SMRT Analysis](http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)  | Self-contained software suite designed for use with Single Molecule, Real-Time (SMRT) Sequencing data.| Analysis Package| [Download](http://www.pacb.com/support/software-downloads)
| [SPAdes](http://bioinf.spbau.ru/en/content/spades-download-0) | Genome assembler intended for both standard isolates and single-cell MDA bacteria assemblies using Illumina or IonTorrent reads and is capable of providing hybrid assemblies using PacBio, Oxford Nanopore and Sanger reads. | Hybrid Assembly | [Download](http://bioinf.spbau.ru/en/content/spades-download-0)
| [Sprai](http://zombie.cb.k.u-tokyo.ac.jp/sprai/README.html)| Sprai (single-pass read accuracy improver) is a tool to correct sequencing errors in single-pass reads for de novo assembly. | Sequencing Error-correction| [Download](http://zombie.cb.k.u-tokyo.ac.jp/sprai/Download.html)


##Illumina Sequencing
###Referenced
| Program        | Description    | Purpose          | Source           |
| :------------- | :------------- | :-------------   | :-------------   |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)| An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. | Reference Mapping | [Download](https://github.com/BenLangmead/bowtie2/releases/latest)
| [BWA](https://github.com/lh3/bwa) | Mapping DNA sequences against a large reference genome, such as the human genome. It consists of three algorithms: **BWA-backtrack**, **BWA-SW** and **BWA-MEM**. | Reference Mapping | [Download](https://github.com/lh3/bwa/releases/latest)

###De novo
| Program        | Description    | Purpose          | Source           |
| :------------- | :------------- | :-------------   | :-------------   |
| [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss) | De novo, parallel, paired-end sequence assembler designed for short reads and large genomes. | Genome Assembly | [Download](https://github.com/bcgsc/abyss/releases/latest) <br> [Install](https://github.com/bcgsc/abyss#quick-start)
| [ALLPATHS-LG](http://www.broadinstitute.org/software/allpaths-lg/blog/) | Short read assembler and it works on both small and large (mammalian size) genomes.| Genome Assembly | [Download](http://www.broadinstitute.org/software/allpaths-lg/blog/?page_id=12)
| [DISCOVAR](http://www.broadinstitute.org/scientific-community/science/programs/genome-sequencing-and-analysis/computational-rd/computational-) | Genome assembler and variant caller. | Genome Assembly | [Download](http://www.broadinstitute.org/software/discovar/blog/?page_id=98)
| [SOAPdenovo](http://soap.genomics.org.cn/soapdenovo.html) | Novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes. | Genome Assembly | [Download](https://sourceforge.net/projects/soapdenovo2/files/latest/download?source=files)
| [SPAdes](http://bioinf.spbau.ru/en/content/spades-download-0) | Genome assembler intended for both standard isolates and single-cell MDA bacteria assemblies using Illumina or IonTorrent reads and is capable of providing hybrid assemblies using PacBio, Oxford Nanopore and Sanger reads. | Genome/Hybrid Assembly | [Download](http://bioinf.spbau.ru/en/content/spades-download-0)
| [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) | Trinity assembles transcript sequences from Illumina RNA-Seq data. | Transcriptome Assembly| [Download](https://github.com/trinityrnaseq/trinityrnaseq/releases/latest)
| [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/) | Short read de novo assembler using de Bruijn graphs. | Genome Assembly | [Download](https://github.com/dzerbino/velvet/tree/master)
