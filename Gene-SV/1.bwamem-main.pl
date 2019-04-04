#!/usr/bin/perl -w
#2017/4/20  

$gene_structure=$ARGV[0];  #"gene.stucture.new";
$ref_genome=$ARGV[1];      # the second genome  (Mo17/B73)
$cds_seq=$ARGV[2]:         #"gene-full-cds-double.new.fa";
$cds_loc=$ARGV[3];         #"cds.loc" 
$script_dir=$ARGV[4];    

`bwa mem -t 1  $ref_genome  gene.fa | samtools view  > aln.sam `; # 'gene.fa' can be splited from 'gene-full-cds-double.new.positive.fa'
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
	#open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
	#`perl $script_dir/2.sam2variation.pl `;
	&sam2variation($prescaffold,$hash{$prescaffold});
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
#open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
#`perl $script_dir/2.sam2variation.pl  `;
&sam2variation($prescaffold,$hash{$prescaffold});
`cat bwasw.out.line >> bwasw.out`;
`rm temp*`;
##################################################
`perl $script_dir/3.check-sequence.pl  $ref_genome    $cds_seq `;
`perl $script_dir/4.clustalw-utl-exon-intron.pl   clustalw.out.checked.filtered  $gene_structure `;
`perl $script_dir/5.Gene_SV_classify.pl    clustalw.out.checked.filtered.structure  `;



#########################################  sam2variation
#########################################
#########################################
sub  sam2variation{

  my ($line1,$line2)=@_;
  my $chr=$line1; my $len=length($line1);  my $temp=" "x(30-$len);  my $scaffoldname="$line1$temp";
  my $ss2=$line2; 

  ####################
  my $find=0;  my $good="";
  my $file="temp";
  open(AA,"$file");
  while(my $line=<AA>)
  { chomp $line;
   my @bb=split/\s+/,$line;
   if(@bb>8){
   my $ll=$bb[5];
   $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
   my @arr=split/\s+/,$ll;
   my $i=0;  my $len=@arr;
   if($arr[1] eq "H" && $arr[0]<50)
     {   if($arr[$len-1] eq "H" && $arr[$len-2]<50){my $find=1;}
		 if($arr[$len-1] eq "M"){my $find=1;}
	 }
   if($arr[$len-1] eq "H" && $arr[$len-2]<50)
    {   if($arr[1] eq "H" && $arr[0]<50){my $find=1;}
	     if($arr[1] eq "M"){my $find=1;}
	}
	if($find>0 && $good eq ""){my $good=$line;}
   }
   my $find=0;
}
close AA;
`rm temp.sam2`;
if($good eq ""){`cp $file temp.sam2`;}
else{open(FA,">temp.sam2");  print FA $good,"\n"; close FA;}
#####################

my $dir0=0; my $dir16=0;
open(AA,"temp.sam2");
open(BB,">temp.out");
while(my $line=<AA>)
{  chomp $line;
   my @bb=split/\s+/,$line;
   if(@bb>8){
   my $ll=$bb[5];
   $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
   my @arr=split/\s+/,$ll;
   my $i=0;  my $len=@arr;
   my $match=0;  my $insertion=0;  my $deletion=0;
   $start1=$bb[3];#
   if($arr[1] eq "H"){my $start2=$arr[0]+1; }else{my $start2=1;}
   while($i<$len)
   {
     if($arr[$i] eq "M"){$match=$match+$arr[$i-1];}
	 if($arr[$i] eq "I"){$insertion=$insertion+$arr[$i-1];}
	 if($arr[$i] eq "D"){$deletion=$deletion+$arr[$i-1];}
	 $i++;
   }
   my $end1=$start1+$match+$deletion-1;
   my $end2=$start2+$match+$insertion-1;
   if($bb[4]>=0) ##not unique
   {
   print BB $bb[0],"\t",$start2,"\t",$end2,"\t",$bb[2],"\t",$start1,"\t",$end1,"\t",$bb[5],"\t",$bb[1],"\n";
  if($bb[1]==0){my $dir0=$dir0+$match;}
    elsif($bb[1]==16){my $dir16=$dir16+$match;}
  }
  }
}
close AA; close BB;
 if($dir0>$dir16){my $direction="+"; my $tt=0;}else{my $direction="-"; my $tt=16;}
open(AA,"temp.out");
open(BB,">temp.out2");
while($line=<AA>)
{
  chomp $line; 
  @bb=split/\s+/,$line;
  if($tt==$bb[7]){print BB "$bb[0]\t$bb[1]\t$bb[2]\t$bb[3]\t$bb[4]\t$bb[5]\t$bb[6]\n";}
}
close AA; close BB;
 
`sort -dk2,2n temp.out2 > temp.out.sorted`;

my $s=0;
open(AA,"temp.out.sorted");
open(BB,">temp.out.sorted.filterred");
while(my $line=<AA>)
{
  chomp $line;
  my @bb=split/\s+/,$line;
  if($s>0){
    if($bb[4]<$start2)
     {if($prelen<($bb[2]-$bb[1])){my $preline=$line; my $prelen=$bb[2]-$bb[1];  my $start1=$bb[1];  my $start2=$bb[4];  my $end1=$bb[2];    my $end2=$bb[5];}}
      elsif($bb[1]<$end1 && $bb[4]<$end2)
	   {if($prelen<($bb[2]-$bb[1])){my $preline=$line; my $prelen=$bb[2]-$bb[1];  my $start1=$bb[1];  my $start2=$bb[4];  my $end1=$bb[2];    my $end2=$bb[5];}}
	  else{ print BB $preline,"\n"; my $preline=$line; my $prelen=$bb[2]-$bb[1];  my $start1=$bb[1];  my $start2=$bb[4];  my $end1=$bb[2];    my $end2=$bb[5];}
   }
   else{my $preline=$line; my $prelen=$bb[2]-$bb[1];  my $start1=$bb[1];  my $start2=$bb[4];  my $end1=$bb[2];    my $end2=$bb[5];}
   $s=1;
}
print BB $preline,"\n";
close AA; close BB;

###################################################################################################################################
if(my $direction eq "+"){my $file="gene1.fa";}else{my $file="gene2.fa";}
open(AA,$file);
my $line1=<AA>;  chomp $line1; $line1=~s/>//;   $line1=~s/\r/ /; my $len=length($line1);  my $temp=" "x(30-$len);  my $genename="$line1$temp";
my $line2=<AA>;  my $ss1=$line2;  chomp $ss1;   my  $genelen=length($ss1); 
close AA;
###################################################################################################################################

my $small=10000000000; my  $big=0;
my $small2=10000000000; my  $big2=0;
open(CC,"temp.out.sorted.filterred");
open(DD,">temp.out.sorted.result");
while(my $line=<CC>)
{
  chomp $line;
  my @bb=split/\s+/,$line;    
  if($bb[1]<$small){$small=$bb[1];}    if($bb[2]>$big){$big=$bb[2];} 
  if($bb[4]<$small2){$small2=$bb[4];}  if($bb[5]>$big2){$big2=$bb[5];} 
  my $len1=$bb[2]-$bb[1]+1;
  my $seq1=substr($ss1,$bb[1]-1,$len1);
  my @str1=split//,$seq1;
  
  my $len2=$bb[5]-$bb[4]+1;
  my $seq2=substr($ss2,$bb[4]-1,$len2);
  my @str2=split//,$seq2;
  
  
  my $ll=$bb[6];
  $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
  my @arr=split/\s+/,$ll;
  my $i=0;  my $len=@arr;   my $sum1=0;   my $sum2=0;  my $new1=""; my $new2=""; 
  while($i<$len)
   {
     if($arr[$i] eq "M")
	   {  my $temp1=substr($seq1,$sum1,$arr[$i-1]);  my $new1="$new1$temp1";   my $sum1=$sum1+$arr[$i-1];  
	      my $temp2=substr($seq2,$sum2,$arr[$i-1]);  my $new2="$new2$temp2";   my $sum2=$sum2+$arr[$i-1];
	   }
	 if($arr[$i] eq "I")
	   {
	      my $temp1=substr($seq1,$sum1,$arr[$i-1]);  my $new1="$new1$temp1";	 my $sum1=$sum1+$arr[$i-1];
	      my $temp2="-"x$arr[$i-1];  my $new2="$new2$temp2";
	   }
	 if($arr[$i] eq "D")
	   { 
	      my $temp1="-"x$arr[$i-1];  my $new1="$new1$temp1";
		  my $temp2=substr($seq2,$sum2,$arr[$i-1]);  my $new2="$new2$temp2";   my $sum2=$sum2+$arr[$i-1];  
	   }
	 
	 $i++;
   }
   print DD $line,"\t",$new1,"\t",$new2,"\n";
}
close CC;  close DD;
######
my $len1=$small-1;
my $addgene1=substr($ss1,0,$len1);  my $addgene11="-"x$len1;
my $len2=$genelen-$big;
my $addgene2=substr($ss1,$big,$len2); my  $addgene22="-"x$len2;
#print $big,"\t",$len2,"\t",$genelen,"\n";


#############
my $sting1=""; my  $string2="";
my $start1=0; my  $start2=0;
my $end1=0; my  $end2=0;
open(CC,"temp.out.sorted.result");
open(DD,">temp.out.sorted.result2");
while(my $line=<CC>)
{
  chomp $line;
  my @bb=split/\s+/,$line;  
  if($bb[1]>$end1 && $bb[4]>$end2)
    {  my $len1=abs($bb[1]-$end1)-1;
	   my $len2=abs($bb[4]-$end2)-1;
	   my $indel1=substr($ss1,$end1,$len1);  $add1="-"x$len1;  
	  my  $indel2=substr($ss2,$end2,$len2);  $add2="-"x$len2;  #print $bb[1],"\t",$end1,"\t",$len1,"\t",$bb[4],"\t",$end2,"\t",$len2,"\n";
	   if($end1==0 && $end2==0)my {$string1="$bb[7]"; my    $string2="$bb[8]";}
	   else{ my $string1="$string1$indel1$add2$bb[7]"; my   $string2="$string2$add1$indel2$bb[8]";}
	}
    elsif($bb[1]<=$end1 && $bb[4]>$end2)
	{ my  $len1=abs($bb[1]-$end1)+1;
	  my  $len2=abs($bb[4]-$end2)-1;
	  my  $indel2=substr($ss2,$end2,$len2);  $add2="-"x$len2;
	  my  $string1="$string1$add2";  $string2="$string2$indel2";
	  my  @temparr=split//,$bb[7]; my  $len_temparr=@temparr;  my  $ii=0; my  $cc="-";
	   while($ii<$len1)
	   {  my  $string1="$string1$cc";  $ii++;}
	   while($ii<$len_temparr)
	   { my  $string1="$string1$temparr[$ii]";  $ii++;}
	   my $string2="$string2$bb[8]";
	}
    elsif($bb[1]>$end1 && $bb[4]<=$end2)
	{ my   $len1=abs($bb[1]-$end1)-1;
	  my  $len2=abs($bb[4]-$end2)+1;
	  my  $indel1=substr($ss1,$end1,$len1);  $add1="-"x$len1;
	  my  $string1="$string1$indel1";
	  my  $string1="$string1$bb[7]";
	  my  @temparr=split//,$bb[8]; my  $len_temparr=@temparr; my  $ii=0;
	  my  $string2="$string2$add1"; $cc="-";
	   while($ii<$len2)
	   { my $string2="$string2$cc";  $ii++;}
	   while($ii<$len_temparr)
	   { my $string2="$string2$temparr[$ii]";  $ii++;}
	}
	my $start1=$bb[1];  my $start2=$bb[4];
    my $end1=$bb[2];   my  $end2=$bb[5];
}
print DD  $addgene1,$string1,$addgene2,"\n",$addgene11,$string2,$addgene22,"\n";
close CC;  close DD;

open(DD,"temp.out.sorted.result2");
my $line1=<DD>; my $line2=<DD>;
close DD;
chomp $line1; chomp $line2;
my $match="";
my @a1=split//,$line1;
my @a2=split//,$line2;
my $len=@a1;
my $i=0;
while($i<$len)
{
 if($a1[$i] eq $a2[$i])
 {
  my $match="$match*";
 }
 else{my $match="$match ";}
 $i++;
}
my $temp=" "x30;
open(EE,">bwasw.out.line");
print EE  ">$scaffoldname\t$genename\t$direction\t$chr\t$small2\t$big2\n";
print EE  "$genename$line1\n";
print EE  "$scaffoldname$line2\n";
print EE  "$temp$match\n";
close EE;

}
#####################################
#####################################
#####################################
