#! /usr/bin/perl -w

# QnD script, Mathias Vandenbogaert, December 2023

#use strict;
use warnings;
use File::Basename;

#recupere les arguments
my ($fastqfile) = @ARGV;

open(FASTQFILE, "<$fastqfile") || die "$fastqfile non trouve\n";

my %all_read_id_list;
my $i=0;
my $one_nb=0;
my $two_nb=0;
my $seq='';
my $read='';

my $run_id=`head -n 1 $fastqfile`;
chomp($run_id);
$run_id =~ s/^@//;
$run_id =~ s/\/[12]$//;
$run_id =~ s/[:\.\s]+.*$//;
print "run_id is $run_id \n";
my %all_read_ids = ();
my $num_orig_seqs = 0;
my $which_pair = 0; # 1 or 2
while (<FASTQFILE>) {
	my $line = $_;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	if ($line =~ m/\@$run_id/) {#illumina-fastq
		if ($which_pair eq "1"){
			$reads_one{$read} = $seq;
			$one_nb++;
		}
		elsif($which_pair eq "2"){
			$reads_two{$read} = $seq;
			$two_nb++;
		}
		if ($read ne ""){
			$all_read_ids{$read} = 1;
		}
		$read = $line;
		($which_pair) = $read =~ m/.*\/([12])$/;
		#print STDERR "which_pair is $which_pair  \n";
		$read =~ s/^@//;
                $read =~ s/\/[12]$//;
		$seq="$line\n";
	        $num_orig_seqs++;
	}
	else {
		$seq.="$line\n";
	}
}
		if ($which_pair eq "1"){
			$reads_one{$read} = $seq;
			$one_nb++;
		}
		elsif($which_pair eq "2"){
			$reads_two{$read} = $seq;
			$two_nb++;
		}

print STDERR "Number of sequences in the original file:\t1:$one_nb -- 2:$two_nb \n";

my $outfile_prefix = basename($fastqfile, ".fastq");
$outfile_prefix = basename($outfile_prefix, ".fq");
# print STDERR "outfile_prefix is $outfile_prefix \n";

open(FONE, ">$outfile_prefix"."_R1.fq") or die;
open(FTWO, ">$outfile_prefix"."_R2.fq") or die;
open(FSGL, ">$outfile_prefix"."_sgl.fq") or die;

my $out_one_nb = 0;
my $out_two_nb = 0;
my $out_sgl_nb = 0;

foreach my $id (sort keys %all_read_ids){
	if (exists($reads_one{$id}) and exists($reads_two{$id})){
	print FONE $reads_one{$id} . "\n";
	print FTWO $reads_two{$id} . "\n";
	$out_one_nb++;
	$out_two_nb++;
}
else{
if (exists($reads_one{$id}) and !exists($reads_two{$id})){
	print FSGL $reads_one{$id} . "\n";
	$out_sgl_nb++;
}
elsif(exists($reads_two{$id}) and !exists($reads_one{$id})){
	print FSGL $reads_two{$id} . "\n";
	$out_sgl_nb++;
}

	
}
}
print STDERR "Number of sequences in the output file:\t\t1:$out_one_nb -- 2:$out_two_nb -- sgl:$out_sgl_nb \n";


close (FASTQFILE);
close(FONE);
close(FTWO);
close(FSGL);
