[![Github Issues](http://githubbadges.herokuapp.com/bretonics/Bioinformatics-Software/issues.svg)](https://github.com/bretonics/Bioinformatics-Software/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/Bioinformatics-Software/pulls.svg)](https://github.com/bretonics/Bioinformatics-Software/pulls)
![](https://reposs.herokuapp.com/?path=bretonics/Bioinformatics-Software&color=orange)


#Tools
| Program        | Description    | Download         |
| :------------- | :------------- | :-------------   |
|[BLAST+] (https://www.ncbi.nlm.nih.gov/books/NBK279690/) | Command line application suite of BLAST tools that utilizes the NCBI C++ Toolkit. | [Download] (ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
| [E-utilities] (http://www.ncbi.nlm.nih.gov/books/NBK25497/) | Entrez Programming Utilities (E-utilities) are a set of nine server-side programs that provide a stable interface into the Entrez query and database system at the NCBI. | N/A
| [EDirect] (http://www.ncbi.nlm.nih.gov/books/NBK179288/) | An advanced method for accessing the NCBI's set of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a UNIX terminal. | N/A
| [SRA Toolkit] (http://www.ncbi.nlm.nih.gov/books/NBK242621/) | The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives. | [Download] (http://ncbi.github.io/sra-tools/)
| [Seqtk](https://github.com/lh3/seqtk) | Fast and lightweight tool for processing sequences in the FASTA or FASTQ format. | [Git Clone](https://github.com/lh3/seqtk.git)


#Analysis
| Program        | Description    | Purpose          | Download         |
| :------------- | :------------- | :-------------   | :-------------   |
| [Galaxy](https://galaxyproject.org/) | Web portal for accessible, reproducible, and transparent computational research. | Analysis package | [Download](https://wiki.galaxyproject.org/Admin/GetGalaxy)
| [Samtools](http://www.htslib.org/) | A suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:<br><br>**Samtools:** Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format.<br><br>**BCFtools:** Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants.<br><br>**HTSlib:** A C library for reading/writing high-throughput sequencing data. | NGS file Toolkit| [Download](http://www.htslib.org/download/)
| [VCFtools](https://vcftools.github.io/index.html) | Package designed for working with complex genetic variation data in the form of VCF files.| VCF Toolkit| [Download](https://vcftools.github.io/downloads.html)
| [MUMmer Package](http://mummer.sourceforge.net/)| Ultra-fast alignment of large-scale DNA and protein sequences. A system for rapidly aligning entire genomes, whether in complete or draft form. <br><br>**MUMmer** is a suffix tree algorithm designed to find maximal exact matches of some minimum length between two input sequences.<br><br>**NUCmer** is a standard DNA sequence alignment. It is a robust pipeline that allows for multiple reference and multiple query sequences to be aligned in a many vs. many fashion.<br><br>**PROmer** is like NUCmer with one exception - all matching and alignment routines are performed on the six frame amino acid translation of the DNA input sequence. | Genome Aligner | [Download](https://sourceforge.net/projects/mummer/files/latest/download?source=files)
| [Mauve](http://darlinglab.org/mauve/mauve.html) | A system for constructing multiple genome alignments in the presence of large-scale evolutionary events such as rearrangement and inversion. | Genome Aligner| [Download](http://darlinglab.org/mauve/download.html)
| [QUAST](http://bioinf.spbau.ru/quast) | Evaluates genome assemblies. | Evaluate Genome Assemblies | [Download](https://sourceforge.net/projects/quast/files/latest/download?source=files)
| [HMMER]() | Search sequence databases for sequence homologs, and for making sequence alignments, analyzed by using profile hidden Markov models | Detect Homologs | [Download](http://hmmer.org/download.html)
| [MUSCLE](http://www.drive5.com/muscle/) | MUSCLE can align hundreds of sequences in seconds. | MSA | [Download]()
| [Clustal](http://www.clustal.org/) | Fast and scalable multiple sequence alignment (can align hundreds of thousands of sequences in hours)| MSA | [Download](http://www.clustal.org/omega/#Download)
| [T-Coffee](http://www.tcoffee.org/) | A multiple sequence alignment package that can align sequences (Protein, DNA, and RNA) or combine the output of your favorite alignment methods (Clustal, Mafft, Probcons, Muscle...) into one unique alignment (**M-Coffee**). It is also able to combine sequence information with protein structural information (**3D-Coffee/Expresso**), profile information (**PSI-Coffee**) or RNA secondary structures. | MSA | [Download](http://www.tcoffee.org/Projects/tcoffee/#DOWNLOAD)
| [QIIME](http://qiime.org/) | Bioinformatics pipeline for performing microbiome analysis from raw DNA sequencing data. QIIME is designed to take users from raw sequencing data generated on the Illumina or other platforms through publication quality graphics and statistics. | Microbial Ecology Pipeline | [Install](http://qiime.org/install/index.html)
| [Mothur](http://www.mothur.org/) |  Expandable software to fill the bioinformatics needs of the microbial ecology community. | Microbial Ecology Pipeline | [Download](http://www.mothur.org/wiki/Download_mothur)

###PacBio Sequencing
| Program        | Description    | Purpose          | Download         |
| :------------- | :------------- | :-------------   | :-------------   |
| [SMRT Analysis](http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)  | Self-contained software suite designed for use with Single Molecule, Real-Time (SMRT) Sequencing data.| Analysis Package| [Download] (http://www.pacb.com/support/software-downloads)
| [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page)| Celera Assembler is a de novo whole-genome shotgun (WGS) DNA sequence assembler, and can use any combination of platform reads. | Genome Assembly | [Download] (https://sourceforge.net/projects/wgs-assembler/files/latest/download?source=files)
| [Cerulean](https://sourceforge.net/projects/ceruleanassembler/) | Cerulean extends contigs assembled using short read datasets like Illumina paired-end reads using long reads like PacBio RS long reads. | Hybrid Assembly| [Download] (https://sourceforge.net/projects/ceruleanassembler/files/latest/download)
| [PBSuite](https://sourceforge.net/projects/pb-jelly/) | **PBJelly** is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes. <br><br>**PBHoney** is an implementation of two variant-identification approaches designed to exploit the high mappability of long reads (i.e., greater than 10,000 bp). PBHoney considers both intra-read discordance and soft-clipped tails of long reads to identify structural variants. | Read Aligner Assembly<br><br>Variant Calling| [Download] (https://sourceforge.net/projects/pb-jelly/files/latest/download)
| [Sprai](http://zombie.cb.k.u-tokyo.ac.jp/sprai/README.html)| Sprai (single-pass read accuracy improver) is a tool to correct sequencing errors in single-pass reads for de novo assembly. | Sequencing error-correction| [Download](http://zombie.cb.k.u-tokyo.ac.jp/sprai/Download.html)


###Illumina Sequencing
####Referenced
| Program        | Description    | Purpose          | Download         |
| :------------- | :------------- | :-------------   | :-------------   |
| [Bowtie](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)| An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. | Read Aligner Assembly| [Download](https://sourceforge.net/projects/bowtie-bio/files/latest/download?source=files)
| [BWA](https://github.com/lh3/bwa) | Mapping DNA sequences against a large reference genome, such as the human genome. It consists of three algorithms: **BWA-backtrack**, **BWA-SW** and **BWA-MEM**. | Read Aligner Assembly | [Download](https://sourceforge.net/projects/bio-bwa/files/latest/download?source=files)

####De novo
| Program        | Description    | Purpose          | Download         |
| :------------- | :------------- | :-------------   | :-------------   |
| [SOAPdenovo](http://soap.genomics.org.cn/soapdenovo.html) | Novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes. | Genome Assembly | [Download] (https://sourceforge.net/projects/soapdenovo2/files/latest/download?source=files)
| [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) | Trinity assembles transcript sequences from Illumina RNA-Seq data. | Transcriptome Assembly| [Download](https://github.com/trinityrnaseq/trinityrnaseq/releases/latest)
| [DISCOVAR](http://www.broadinstitute.org/scientific-community/science/programs/genome-sequencing-and-analysis/computational-rd/computational-) | Genome assembler and variant caller. | Genome Assembly | [Download](ftp.broadinstitute.org/pub/crd/Discovar/)
| [ALLPATHS-LG](http://www.broadinstitute.org/software/allpaths-lg/blog/) | Short read assembler and it works on both small and large (mammalian size) genomes.| Genome Assembly | [Download](ftp.broadinstitute.org/pub/crd/ALLPATHS/Release-LG/latest_source_code/LATEST_VERSION.tar.gz)
| [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/) | Short read de novo assembler using de Bruijn graphs. | Genome Assembly | [Download](https://github.com/dzerbino/velvet/tree/master)
