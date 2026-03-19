#! /usr/bin/perl -w

# QnD script, Mathias Vandenbogaert, December 2023

#programme qui prend en entree
#1 - le fichier contenant la liste des sequences
#2 - le fichier fastq de toutes les séquences
#3 - le fichier de sortie

#le programme cree un fichier fastq avec uniquement les sequences d'interet

#use strict;
use warnings;

#recupere les arguments
my ($infile, $fastqfile, $outfile) = @ARGV;

open(INFILE, "<$infile") || die "$infile non trouve\n";
open(FASTQFILE, "<$fastqfile") || die "$fastqfile non trouve\n";
open(OUTFILE, ">$outfile");

my %list;
my $i=0;
my $seq='';
my $contig='';

my $run_id=`head -n 1 $fastqfile`;
chomp($run_id);
$run_id =~ s/^@//;
$run_id =~ s/\/[12]$//;
$run_id =~ s/[:\.\s]+.*$//;
print "run_id is $run_id \n";

my $num_ids = 0;
while (<INFILE>) {
	my $line = $_;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	my @param = split (/\t/,$line);
	$list{$param[0]}=1;
	$num_ids++;
}
close (INFILE);
print STDERR "Number of query IDs: $num_ids \n";

my $num_orig_seqs = 0;
while (<FASTQFILE>) {
	my $line = $_;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	if ($line =~ m/\@$run_id/) {#illumina-fastq
		if (exists $list{$contig}) {
			print OUTFILE "$seq";
			$i++;
		}
		$contig = $line;
		$contig =~ s/^@//;
                $contig =~ s/\/[12]$//;
		$seq="$line\n";
	        $num_orig_seqs++;
	}
	else {
		$seq.="$line\n";
	}
}

if (exists $list{$contig}) {
	print OUTFILE "$seq";
	$i++;
}
print STDERR "Number of sequences in the original file: $num_orig_seqs \n";
print STDERR "Number of output sequences: $i \(= " . $i/2 . " * 2\)\n";

close (FASTQFILE);
close (OUTFILE);
