#!/usr/bin/perl -w
#2017/4/18  yszhou

open(AA,"/share/home/caulai/yszhou/data/B73-genome-v4/merged-annotation/B73-merged-filtered.gff3");
open(BB,">B73-cds.loc");
while($line=<AA>)
{ chomp $line;
  $line=~s/ID=CDS://;  $line=~s/;/\t/; 
  @bb=split/\s+/,$line;
  if($bb[2] eq "CDS")
  {
    if($pre ne $bb[8])
	  {
	    if($s>0){print BB "$chr\t$start\t$end\t$dir\t$pre\n";}
		$chr=$bb[0]; $start=$bb[3]; $end=$bb[4]; $dir=$bb[6]; $pre=$bb[8];
	  }
	  else{
	    if($bb[3]<$start){$start=$bb[3];}
		if($bb[4]>$end){$end=$bb[4];}
	  }
	  $s=1;
  }
}
print BB "$chr\t$start\t$end\t$dir\t$pre\n";
close AA; close BB;
###################################################################################################
%hash=();  %lencht=();
open(FF,"/share/home/caulai/yszhou/data/B73-genome-v4/maize_pseudo4-full.fa");
while($line=<FF>)
{ chomp $line;
  if($line=~/^>/)
  { @bb=split/\s+/,$line;  $chr=$bb[0]; $chr=~s/>//;}
  else{$hash{$chr}=$line;
  $lenchr{$chr}=length($line);
  }
}
close FF;

open(CC,"B73-cds.loc");
open(DD,">B73-gene-full-cds-2000.fa");
while($line=<CC>)
{ chomp $line;
  @bb=split/\s+/,$line;
  if($bb[1]>2000){$nu1=2000;}else{$nu1=$bb[1]-1;}
  if(($bb[2]+2000)<$lenchr{$bb[0]}){$nu2=2000;}else{$nu2=$lenchr{$bb[0]}-$bb[2]-1;}
  
  $len=abs($bb[2]-$bb[1])+1+$nu1+$nu2;
  $seq=substr($hash{$bb[0]},$bb[1]-1-$nu1,$len);
  if($bb[3] eq "+")
  {print DD ">$bb[4]\t$nu1\t$nu2\n$seq\n"; }
  else{
   $len=length($seq); @ar=split//,$seq;
   $i=0;  $seq2="";
	while($i<$len)
	{$seq2="$ar[$i]$seq2";
	 $i++;
	}
	$seq2=~s/A/x/g;  $seq2=~s/T/y/g;  $seq2=~s/x/T/g;  $seq2=~s/y/A/g;  
	$seq2=~s/C/x/g;  $seq2=~s/G/y/g;  $seq2=~s/x/G/g;  $seq2=~s/y/C/g;  
	 print DD ">$bb[4]\t$nu2\t$nu1\n$seq\n";
  }
}
close CC; close DD;

##################################################################################################
%start=();  %end=();  %dir=();
$file1="/share/home/caulai/yszhou/data/B73-genome-v4/merged-annotation/B73-merged-filtered.gff3";  ##input: gff3 
$file2="B73.gene.stucture";  ## output: gene simple  structure
open(FF,"B73-cds.loc");
while($line=<FF>)
{ chomp $line;
  @bb=split/\s+/,$line;
  $start{$bb[4]}=$bb[1];
  $end{$bb[4]}=$bb[2];
  $dir{$bb[4]}=$bb[3];
}
close FF;
##############################
$s=0;
open(AA,"$file1");
open(BB,">$file2");
while($line=<AA>)
{  chomp $line;
   $line=~s/Parent=//;  $line=~s/ID=CDS\://;  $line=~s/ID=//;  $line=~s/;/\t/;  
   @bb=split/\s+/,$line;
   ###########################
   if($bb[2] eq "CDS")
   {  if($pre ne $bb[8])
       {  if($s>0){print BB ">$pre\t$dir{$pre}\n",$total,"\n";}
          if($bb[6] eq "+")
		  { $t1=$bb[3]-$start{$bb[8]}+1; $t2=$bb[4]-$start{$bb[8]}+1;}
		  else{$t1=$end{$bb[8]}-$bb[4]+1;  $t2=$end{$bb[8]}-$bb[3]+1;}
		  $total="EXON\t$t1\t$t2";
		  $s++;
	   }
	   else{   if($bb[6] eq "+"){
	                $t3=$bb[3]-$start{$bb[8]}+1; $t4=$bb[4]-$start{$bb[8]}+1;
	                $intron1=$t2+1; $intron2=$t3-1;
					$total="$total\nINTRON\t$intron1\t$intron2\nEXON\t$t3\t$t4";
					}
				else{
				   $t3=$end{$bb[8]}-$bb[4]+1;  $t4=$end{$bb[8]}-$bb[3]+1;
				   $intron1=$t4+1; $intron2=$t1-1;
				   $total="EXON\t$t3\t$t4\nINTRON\t$intron1\t$intron2\n$total";
				}
					$t1=$t3; $t2=$t4;	
	       }
		$pre=$bb[8]; 
	}
}
print BB ">$pre\t$dir{$pre}\n",$total,"\n";
close AA;  close BB;
#######################################################################################
%h1=();  %h2=();
open(AA,"B73-gene-full-cds-2000.fa");
while($line=<AA>)
{ chomp $line;
  if($line=~/^>/)
  {
   $line=~s/>//;
   @bb=split/\s+/,$line;
   $n1{$bb[0]}=$bb[1];
   $n2{$bb[0]}=$bb[2];
  }
}
close AA;

$s=0;
open(BB,"B73.gene.stucture");
open(CC,">B73.gene.stucture.new");
while($line=<BB>)
{ chomp $line;
  if($line=~/^>/)
  {  if($s>0){print CC "UTL2\t$pre1\t$pre2\n";}
      $ll=$line;  $ll=~s/>//;
      @bb=split/\s+/,$ll;
      $nu1=$n1{$bb[0]};
	  $nu2=$n2{$bb[0]};
	  print CC $line,"\n";
	  print CC "UTL1\t","1\t$nu1\n";
    $s=1;
  }
  else{  @bb=split/\s+/,$line;
         $t1=$bb[1]+$nu1;
         $t2=$bb[2]+$nu1;
         print CC "$bb[0]\t$t1\t$t2\n";
         $pre1=$bb[2]+$nu1+1; $pre2=$bb[2]+$nu1+$nu2;		 
  }
}
print CC "UTL2\t$pre1\t$pre2\n";
close BB; close CC;
#######################################################################################
%ge=(); %ge2=();
open(AA,"B73-gene-full-cds-2000.fa");
open(BB,">B73-gene-full-cds-double.fa");
while($line=<AA>)
{ chomp $line;
  if($line=~/^>/)
  {$ll=$line; $ll=~s/\|/\t/; @bb=split/\s+/,$ll; $pre=$bb[0]; $pre=~s/>//;}
  else{
   $ge{$pre}=$line;     print BB ">$pre\t1\n$line\n";
     $ll=$line; 
	 $ll=~s/A/x/g;  $ll=~s/T/y/g;  $ll=~s/C/z/g;  $ll=~s/G/k/g; 
	 $ll=~s/x/T/g;  $ll=~s/y/A/g;  $ll=~s/z/G/g;  $ll=~s/k/C/g; 
     @ar=split//,$ll;  $len=@ar;
	 $i=0; $ll2="";
	 while($i<$len)
	 {$ll2="$ar[$i]$ll2"; $i++;}
      $ge2{$pre}=$ll2;
     print BB ">$pre\t2\n$ll2\n";
	 }
}
close AA; close BB;
################################################################################
%hash1=();  %hash2=();
open(AA,"B73-gene-full-cds-double.fa");
while($line=<AA>)
{ chomp $line;
  if($line=~/^>/){$line=~s/>//; @bb=split/\s+/,$line;}
  else{
    if($bb[1]==1){$hash1{$bb[0]}=$line;}
	else{$hash2{$bb[0]}=$line;}
  }
}
close AA;

open(BB,"B73-cds.loc");
open(CC,">B73-gene-full-cds-double.new.fa");
open(DD,">B73-gene-full-cds-double.new.positive.fa");
while($line=<BB>)
{
  chomp $line;
  @bb=split/\s+/,$line;
  if(exists($hash1{$bb[4]})){
  if($bb[3] eq "+")
  {  print CC ">$bb[4]\t1\n",$hash1{$bb[4]},"\n";
     print CC ">$bb[4]\t2\n",$hash2{$bb[4]},"\n";
     print DD ">$bb[4]\t1\n",$hash1{$bb[4]},"\n";
  }
  else{
     print CC ">$bb[4]\t1\n",$hash2{$bb[4]},"\n";
     print CC ">$bb[4]\t2\n",$hash1{$bb[4]},"\n";
     print DD ">$bb[4]\t1\n",$hash2{$bb[4]},"\n";
  }
  }
}
close BB; close CC; close DD;
