#!/usr/bin/perl
use warnings;
use strict;
die "USAGE:perl $0 genome.fa.ltrdigest_5ltr.fas genome.fa.ltrdigest_3ltr.fas genome.fa.ltrdigest_tabout.csv\n" if(@ARGV!=3);

my %lutr;
open A,$ARGV[0] or die $!;
while(<A>)
{
    chomp;
    my $id;
    if(/>/)
    {
        $id=$_;
    }
    else
    {
        $lutr{$id}.=$_;
    }
}
close A;

my %rutr;
open B,$ARGV[1] or die $!;
while(<B>)
{
    chomp;
    my $id;
    if(/>/)
    {
        $id=$_;
    }
    else
    {
        $rutr{$id}.=$_;
    }
}
close B;

open C,$ARGV[2] or die $!;
while(<C>)
{
    chomp;
    my @a=split(/\t/);
    next if($a[29]!~/\w/);
    my $id=">".$a[3]."_".$a[0]."_".$a[1];
    if(exists $lutr{$id})
    {
        open TMP1,">mo17_temp_lutr.fa" or die $!;
        my $lutr_name=">temp_lltr";
        print TMP1 "$lutr_name\n$lutr{$id}\n";
        close TMP1;
    }
    if(exists $rutr{$id})
    {
        open TMP2,">mo17_temp_rutr.fa" or die $!;
        my $rutr_name=">temp_rltr";
        print TMP2 "$rutr_name\n$rutr{$id}\n";
        close TMP2;
    }
    system("cat mo17_temp_lutr.fa mo17_temp_rutr.fa >mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa");
    system("/share/home/caulai/software/muscle3.8.31_i86linux64 -in mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa -out mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa.aln -clwstrict");
    system("/share/home/caulai/software/EMBOSS-6.6.0/bin/distmat -sequence mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa.aln -nucmethod 2 -outfile mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa.aln.dist");

    my $n=0;
    open DIST,"mo17_temp_lutr.fa_vs_mo17_temp_rutr.fa.aln.dist" or die $!;
    open TIME,">>$ARGV[2].time" or die $!;
    while(<DIST>)
    {
        chomp;
        $n++;
        if($n==9)
        {
            my @a=split;
            my $distance=$a[1];
            if($distance eq "-0.00")
            {
                $distance=0;
            }
            my $time=$distance/2*0.000000013;
            $time=$time/100000000;
            print TIME "$id\t$distance\t$time\n";
        }
    }
    close DIST;
    close TIME;
}
close C;
