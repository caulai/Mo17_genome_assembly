#!/usr/bin/perl
use strict;
use warnings;

die usage() if @ARGV == 0;
my ($bam1,$bam2,$coverage_cutoff,$gap_proportion,$output) = @ARGV;

## the self alignment file
my %hash_unaligned_sequence_self;
my %hash_gap_sequence;
open NEW,"samtools view $bam2 |" or die;
while(<NEW>){
	chomp;
	my @array = split /\s+/;
	my $N_number = $array[9] =~ s/N/N/g;
	my $seq_length = length($array[9]);
	if(($N_number/$seq_length) >= $gap_proportion){
		$hash_gap_sequence{$array[0]} = 1;
	}
	elsif($array[1] == 4){
		$hash_unaligned_sequence_self{$array[0]} = 1;
	}
}
close NEW;

## open the cross alignment file
open NEW,"samtools view $bam1 |" or die;
## the output file pf bed file
open NEW1,">$output" or die;
while(<NEW>){
	chomp;
	my @array = split /\s+/;
	$array[0] =~ /(\w+):(\d+)-(\d+)/;
	my $scaffold_name = $1;
	my $start = $2;
	my $end = $3;
	if(exists $hash_gap_sequence{$array[0]} or exists $hash_unaligned_sequence_self{$array[0]}){
		next;
	}
	elsif($array[1] == 4){
		print NEW1 "$scaffold_name\t$start\t$end\n";
	}
	elsif($array[1] == 0 or $array[1] == 16){
		my $match_bases = 0;
		while($array[5] =~ /(\d+)M/g){
			$match_bases+=$1;
		}
		my $seq_length = $end - $start + 1;
		my $coverage = $match_bases/$seq_length;
		if($coverage < $coverage_cutoff){
			print NEW1 "$scaffold_name\t$start\t$end\n";
		}
	}
}

close NEW;
close NEW1;
	
sub usage{
	my $die =<<DIE;
	usage : perl *.pl Mo17_hybrid_final.500_100.bwa_mem.Zea_mays.AGPv4.dna.toplevel.bam Mo17_hybrid_final.500_100.bwa_mem.Mo17_hybrid_final.bam coverage_cutoff gap_proportion output.bed
DIE
}
