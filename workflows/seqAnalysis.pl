#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

# =============================================================================
#
#   CAPITAN: Andres Breton, http://andresbreton.com
#   FILE: seqAnalysis.pl
#   USAGE: Automated Sequence Analysis Workflow example
#   DEPENDENCIES: Software outlined in Sequence Analysis Workflow
#                 (https://github.com/bretonics/Bioinformatics-Software/wiki/Sequencing-Analysis-Workflow#software)
#
# =============================================================================

#-------------------------------------------------------------------------
# COMMAND LINE
my $SRA = "";
my $ACCESSION = "";
my $usage = "\n\n$0 [options]\n
Options:
    -sra          SRA ID
    -accession    Reference accession number
    -help         Show this message
\n";

# OPTIONS
GetOptions(
    'sra:s'         =>\$SRA,
    'accession:s'   =>\$ACCESSION,
    help            =>sub{pod2usage($usage);}
)or pod2usage(2);
#-------------------------------------------------------------------------
# VARIABLES
my $index = "MG1655"; #reference index name
my $sraRepo = "~/ncbi/public/sra"; #default SRA directory
my $outDir 	= setOutputDir('analysis'); #create analysis directory

# Hash of Hashes of software with respective commands used in the workflow.
# This simplifies both visual and future modifications you'd like to
# make for each command. One place, one change.
my %COMMANDS = (
    edirect             => "esearch -db nucleotide -query \"$ACCESSION [ACCN]\" | efetch -format fasta > $ACCESSION.fasta",
    prefetch            => "prefetch -c $SRA",
    "fastq-dump"        => "fastq-dump --split-files $SRA ; fastq-dump --fasta $SRA",
    "bowtie2-build"     => "bowtie2-build -o 3 $ACCESSION.fasta $index",
    bowtie2             => "bowtie2 -x $index -1 $SRA\_1.fastq -2 $SRA\_2.fastq -S $SRA.sam",
    samtools            => {
                            view        => "samtools view -b $SRA.sam > $SRA.bam",
                            sort        => "samtools sort $SRA.bam $SRA.sorted", #server has different version, thus no < -o > flag
                            index       => "samtools index $SRA.sorted.bam",
                            tview       => "samtools tview $SRA.sorted.bam",
                            consensus   => "samtools mpileup -uf $ACCESSION.fasta $SRA.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > consensus.fastq",
                            variants    => "samtools mpileup -uf $ACCESSION.fasta $SRA.sorted.bam | bcftools call -mv -Oz > variants.vcf.gz",
                            },
    seqtk               => "seqtk seq -A consensus.fastq > consensus.fasta",
    nucmer              => "nucmer -maxmatch -c 100 -p nucmer $ACCESSION.fasta consensus.fasta",
    "show-coords"       => "show-coords -r -c -l nucmer.delta > nucmer.coords",
    mummerplot          => "mummerplot --png --color -p nucmer nucmer.delta",
);
#-------------------------------------------------------------------------
# CALLS
argChecks(); #check CL arguments were passed
executeCommand("edirect"); #edirect
executeCommand("prefetch"); #prefetch
executeCommand("fastq-dump"); #fastq-dump
executeCommand("bowtie2-build"); #bowtie2-build
executeCommand("bowtie2"); #bowtie2
executeCommand("samtools", "view"); #samtools view
executeCommand("samtools", "sort"); #samtools sort
executeCommand("samtools", "index"); #samtools index
# samtools 'tview' commented b/c you don't want to interupt your analysis
# with a screen...execute command above once done to see alignment
# executeCommand("samtools", "tview"); #samtools tview
executeCommand("samtools", "consensus"); #samtools consensus
executeCommand("seqtk"); #seqtk
executeCommand("nucmer"); #nucmer
executeCommand("show-coords"); #show-coords
executeCommand("mummerplot"); #mummerplot
executeCommand("samtools", "variants"); #samtools variants
#-------------------------------------------------------------------------
# SUBS

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = argChecks();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function checks command-line arguments using global variables
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = Dies from errors
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub argChecks {
    unless ($SRA) {
        die "You did not provide an SRA ID.", $!, $usage;
    }
    unless ($ACCESSION) {
        die "You did not provide an accession number for your reference.", $!, $usage;
    }
    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = executeCommand($call, $option);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes arguments 1 or 2 arguments; the call to
# the program to be executed and an optional option for commands
# requiring a second flag (such as samtools)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = Executes command and reports status
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub executeCommand {
    my ($call, $option) = @_;
    if ($option) {
        say "Executing $call...\n$COMMANDS{$call}{$option}\n";
        my $result = `$COMMANDS{$call}{$option}`;
        failedEx($COMMANDS{$call}{$option}) if ($? != 0);
        say "Done.\n\n";
    } else {
        say "Executing $call...\n$COMMANDS{$call}\n";
        my $result = `$COMMANDS{$call}`;
        failedEx($COMMANDS{$call}) if ($? != 0);
        say "Done.\n\n";
    }
    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($command);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 1 argument, the executed sys call command.
# It prints warnings and dies when command execution fails
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = Print warnings and die when command execution fails
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub failedEx {
    my ($command) = @_;
    die	"WARNING: Something seems to have gone wrong!\n",
    	"Failed to execute command " , $command , "\n\n",
    	"Please check installed software version and/or permissions", $!;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = VOID context
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes no arguments, creates an output directory
# if non-existent for results, and moves into it
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = 'analysis' output directory
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub setOutputDir {
    my ($outDir) =  @_;
    if (! -e $outDir){
        `mkdir $outDir`;
    }
    say "Changing to $outDir directory...\n";
    chdir($outDir) or die "$!";
    return $outDir;
}
