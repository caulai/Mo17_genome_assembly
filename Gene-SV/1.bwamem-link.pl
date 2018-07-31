#!/usr/bin/perl -w
#2017/4/20  

$gene_structure=$ARGV[0];  #"gene.stucture.new";
$ref_genome=$ARGV[1];      # the second genome  (Mo17/B73)
$cds_seq=$ARGV[2]:         #"gene-full-cds-double.new.fa";
$cds_loc=$ARGV[3];         #"cds.loc" 
$script_dir=$ARGV[4];      

`bwa mem -t 1  $ref_genome   gene.fa > aln.sam `;
###################
open(AA,"$ref_genome");
while($line1=<AA>)
{ $line2=<AA>;  
  chomp $line1;  
  chomp $line2;
  $line1=~s/>//;
  @bb=split/\s+/,$line1;  
  $hash{$bb[0]}=$line2;
}close AA;

open(AA,"gene.fa");
while($line=<AA>)
{ chomp $line;
  if($line=~/^>/)
  { $ll=$line; $ll=~s/>//;  @bb=split/\s+/,$ll; $pre=$bb[0];}
  else{
   $ge{$pre}=$line;
   $line2=$line;
   $line2=~tr/ATCGatcg/TAGCtagc/;
   @arr=split//,$line2;
   $lenarr=@arr;
   $i=0;  $st="";
   while($i<$lenarr){$st="$arr[$i]$st"; $i++;}
   $ge2{$pre}=$st;
    }
}
close AA;
##########################################
open(FF,"$cds_loc");
while($line=<FF>)
{  chomp $line;
   @bb=split/\s+/,$line;
   $loc{$bb[4]}=($bb[1]+$bb[2])/2;
   $genechr{$bb[4]}=$bb[0];
}
close FF;

open(AAA,"aln.sam");
open(BBB,">aln.sam.filtered");
while($line=<AAA>)
{
  chomp $line;
  @bb=split/\s+/,$line;
  if($bb[1]==0 || $bb[1]==2048){$dir=0;}
  else{$dir=16;}
  $chr1=$bb[2];    $chr2=$genechr{$bb[0]};
  $chr1=~s/chr//;  $chr2=~s/chr//;
  $imfor=$bb[5];  
  $imfor=~s/S/H/g;
  $lenbb=@bb;
  if($bb[4]>=0 && $lenbb>4)  #not unique
  {  $distance=abs($loc{$bb[0]}-$bb[3]);
    #if($distance<10000000 &&  $chr1 eq $chr2)  #syteny
    {print BBB $bb[0],"\t",$dir,"\t",$bb[2],"\t",$bb[3],"\t",$bb[4],"\t",$imfor,"\t",$bb[6],"\t",$bb[7],"\t",$bb[8],"\t",$bb[9],"\t",$bb[10],"\t",$bb[11],"\t",$bb[12],"\t",$bb[13],"\t",$bb[14],"\n";}
  }
}
close AAA; close BBB;
###########################################
`sort -dk1,1 -k3,3 -k4,4n  aln.sam.filtered  > aln.sam.sorted`;
$s=0;
open(AA,"aln.sam.sorted");
open(FF,">bwasw.out"); close FF;
while($line=<AA>)
{
  chomp $line;
  @bb=split/\s+/,$line;
  if($s==0){ open(BB,">temp");  print BB $line,"\n";  $start=$bb[3];}
  else{
  if($bb[0] eq $pregene && $bb[2] eq $prescaffold  && abs($bb[3]-$prestart)<20000){print BB $line,"\n";}
  else{close  BB;
  open(A1,">gene1.fa");  print A1 ">$pregene\n",$ge{$pregene},"\n";   close A1;
	open(A2,">gene2.fa");  print A2 ">$pregene\n",$ge2{$pregene},"\n";  close A2;
	open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
	`perl $script_dir/2.sam2variation.pl `;
	`cat bwasw.out.line >> bwasw.out`;
	`rm temp*`;
    open(BB,">temp");  print BB $line,"\n";  $start=$bb[3];
  }
  }
  $pregene=$bb[0];  $prescaffold=$bb[2];  $prestart=$bb[3]; $s=1;
}
close AA; 
close BB;
open(A1,">gene1.fa");  print A1 ">$pregene\n",$ge{$pregene},"\n";   close A1;
open(A2,">gene2.fa");  print A2 ">$pregene\n",$ge2{$pregene},"\n";  close A2;
open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
`perl $script_dir/2.sam2variation.pl  `;
`cat bwasw.out.line >> bwasw.out`;
`rm temp*`;
##################################################
`perl $script_dir/3.check-sequence.pl  $ref_genome    $cds_seq `;
`perl $script_dir/4.clustalw-utl-exon-intron.pl  $gene_structure `;
