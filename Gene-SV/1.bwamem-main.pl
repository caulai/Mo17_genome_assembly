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
   #$loc{$bb[4]}=($bb[1]+$bb[2])/2;
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
  {  #$distance=abs($loc{$bb[0]}-$bb[3]);
    #if($distance<10000000 &&  $chr1 eq $chr2)  #syteny
    {print BBB $bb[0],"\t",$dir,"\t",$bb[2],"\t",$bb[3],"\t",$bb[4],"\t",$imfor,"\t",$bb[6],"\t",$bb[7],"\t",$bb[8],"\t",$bb[9],"\t",$bb[10],"\t",$bb[11],"\t",$bb[12],"\t",$bb[13],"\t",$bb[14],"\n";}
  }
}
close AAA; close BBB;
###########################################
`sort -dk1,1 -k3,3 -k4,4n  aln.sam.filtered  > aln.sam.sorted`;
$sr=0;
open(A,"aln.sam.sorted");
open(F,">bwasw.out"); close F;
while($Lines=<A>)
{ 
  chomp $Lines;
  @br=split/\s+/,$Lines;
  if($sr==0){ open(B,">temp");  print B $Lines,"\n";  $startr=$br[3];}
  else{
  if($br[0] eq $pregener && $br[2] eq $prescaffoldr  && abs($br[3]-$prestartr)<20000){print B $Lines,"\n";}
  else{  close B;
  open(A1,">gene1.fa");  print A1 ">$pregener\n",$ge{$pregener},"\n";   close A1;
  open(A2,">gene2.fa");  print A2 ">$pregener\n",$ge2{$pregener},"\n";  close A2;
  #open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
  #`perl $script_dir/2.sam2variation.pl  `;
   &sam2variation($prescaffoldr,$hash{$prescaffoldr});
  `cat bwasw.out.line >> bwasw.out`;
   `rm temp*`;
    open(B,">temp");  print B $Lines,"\n";  $startr=$br[3];
  }
  }
  $pregener=$br[0];  $prescaffoldr=$br[2];  $prestartr=$br[3];   $sr=1;
}
close A; 
close B;
open(A1,">gene1.fa");  print A1 ">$pregene\n",$ge{$pregene},"\n";   close A1;
open(A2,">gene2.fa");  print A2 ">$pregene\n",$ge2{$pregene},"\n";  close A2;
#open(A3,">scaffold.fa");  print A3 ">$prescaffold\n",$hash{$prescaffold},"\n";  close A3;
#`perl $script_dir/2.sam2variation.pl  `;
&sam2variation($prescaffoldr,$hash{$prescaffoldr});
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
   
  my ($line1,$line2,$chr,$len,$temp,$scaffoldname,$s,$ss1,$ss2,$find,$good,$file,$line,$ll,$i,$dir,$dir16,$match,$insertion,$deletion,$start1,$start2,$end1,$end2,$dir0,,$drection,$tt,$preline,$prelen,$seq1,$seq2,$small,$small2,$big,$big2);
  my ($string1,$string2,$indel1,$indel2,$temp1,$temp2,$new1,$new2,$addgene1,$addgene2,$sum1,$sum2,$addgene11,$addgene22,$len_temparr,$ii,$cc,$genename);
  my (@bb,@arr,@str1,@str2,@temparr,@a1,@a2);
   ($line1,$line2)=@_;
   $chr=$line1;  $len=length($line1);   $temp=" "x(30-$len);   $scaffoldname="$line1$temp";
   $ss2=$line2; 

###################
$find=0;  $good="";
$file="temp";
open(AAA,"$file");
while($line=<AAA>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if(@bb>8)
    {
       $ll=$bb[5];
       $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
       @arr=split/\s+/,$ll;
       $i=0;  $len=@arr;
       if($arr[1] eq "H" && $arr[0]<50)
           {   if($arr[$len-1] eq "H" && $arr[$len-2]<50){$find=1;}
	        if($arr[$len-1] eq "M"){$find=1;}
           }
       if($arr[$len-1] eq "H" && $arr[$len-2]<50)
           {   if($arr[1] eq "H" && $arr[0]<50){$find=1;}
	        if($arr[1] eq "M"){$find=1;}
           }
       if($find>0 && $good eq ""){$good=$line;}
   }
   $find=0;
}
close AAA;
#`rm temp.sam2`;
if($good eq ""){`cp $file temp.sam2`;}
else{open(FFA,">temp.sam2");  print FFA $good,"\n";   close FFA;}
#####################

$dir0=0; $dir16=0;
open(AAA,"temp.sam2");
open(BBB,">temp.out");
while($line=<AAA>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if(@bb>8){
   $ll=$bb[5];
   $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
   @arr=split/\s+/,$ll;
   $i=0;  $len=@arr;
   $match=0;  $insertion=0;  $deletion=0;
   $start1=$bb[3];#
   if($arr[1] eq "H"){$start2=$arr[0]+1; }else{$start2=1;}
   while($i<$len)
   {
     if($arr[$i] eq "M"){$match=$match+$arr[$i-1];}
	 if($arr[$i] eq "I"){$insertion=$insertion+$arr[$i-1];}
	 if($arr[$i] eq "D"){$deletion=$deletion+$arr[$i-1];}
	 $i++;
   }
   $end1=$start1+$match+$deletion-1;
   $end2=$start2+$match+$insertion-1;
   if($bb[4]>=0) ##not unique
   {
   print BBB $bb[0],"\t",$start2,"\t",$end2,"\t",$bb[2],"\t",$start1,"\t",$end1,"\t",$bb[5],"\t",$bb[1],"\n";
  if($bb[1]==0){$dir0=$dir0+$match;}
    elsif($bb[1]==16){$dir16=$dir16+$match;}
  }
  }
}
close AAA; close BBB;

if($dir0>$dir16){$direction="+"; $tt=0;}else{$direction="-"; $tt=16;}

open(AAA,"temp.out");
open(BBB,">temp.out2");
while($line=<AAA>)
{
  chomp $line; 
  @bb=split/\s+/,$line;
  if($tt==$bb[7]){print BBB "$bb[0]\t$bb[1]\t$bb[2]\t$bb[3]\t$bb[4]\t$bb[5]\t$bb[6]\n";}
}
close AAA; close BBB;
 
`sort -dk2,2n temp.out2 > temp.out.sorted`;

$s=0;
open(AAA,"temp.out.sorted");
open(BBB,">temp.out.sorted.filterred");
while($line=<AAA>)
{
  chomp $line;
  @bb=split/\s+/,$line;
  if($s>0){
    if($bb[4]<$start2)
     {if($prelen<($bb[2]-$bb[1])){$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}}
      elsif($bb[1]<$end1 && $bb[4]<$end2)
	   {if($prelen<($bb[2]-$bb[1])){$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}}
	  else{ print BBB $preline,"\n"; $preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}
   }
   else{$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}
   $s=1;
}
print BBB $preline,"\n";
close AAA; close BBB;

###################################################################################################################################
if($direction eq "+"){$file="gene1.fa";}else{$file="gene2.fa";}
open(AAA,$file);
$line1=<AAA>;  chomp $line1; $line1=~s/>//;  $line1=~s/\r/ /; $len=length($line1);  $temp=" "x(30-$len);  $genename="$line1$temp";
$line2=<AAA>;  $ss1=$line2;  chomp $ss1;    $genelen=length($ss1); 
close AAA;
###################################################################################################################################

$small=10000000000;  $big=0;
$small2=10000000000;  $big2=0;
open(CCC,"temp.out.sorted.filterred");
open(DDD,">temp.out.sorted.result");
while($line=<CCC>)
{
  chomp $line;
  @bb=split/\s+/,$line;    
  if($bb[1]<$small){$small=$bb[1];}    if($bb[2]>$big){$big=$bb[2];} 
  if($bb[4]<$small2){$small2=$bb[4];}  if($bb[5]>$big2){$big2=$bb[5];} 
  $len1=$bb[2]-$bb[1]+1;
  $seq1=substr($ss1,$bb[1]-1,$len1);
  @str1=split//,$seq1;
  
  $len2=$bb[5]-$bb[4]+1;
  $seq2=substr($ss2,$bb[4]-1,$len2);
  @str2=split//,$seq2;
  
  
  $ll=$bb[6];
  $ll=~s/H/\tH\t/g;  $ll=~s/M/\tM\t/g;  $ll=~s/D/\tD\t/g;   $ll=~s/I/\tI\t/g; 
  @arr=split/\s+/,$ll;
  $i=0;  $len=@arr;   $sum1=0;   $sum2=0;  $new1=""; $new2=""; 
  while($i<$len)
   {
     if($arr[$i] eq "M")
	   {  $temp1=substr($seq1,$sum1,$arr[$i-1]);  $new1="$new1$temp1";   $sum1=$sum1+$arr[$i-1];  
	      $temp2=substr($seq2,$sum2,$arr[$i-1]);  $new2="$new2$temp2";   $sum2=$sum2+$arr[$i-1];
	   }
	 if($arr[$i] eq "I")
	   {
	      $temp1=substr($seq1,$sum1,$arr[$i-1]);  $new1="$new1$temp1";	 $sum1=$sum1+$arr[$i-1];
	      $temp2="-"x$arr[$i-1];  $new2="$new2$temp2";
	   }
	 if($arr[$i] eq "D")
	   { 
	      $temp1="-"x$arr[$i-1];  $new1="$new1$temp1";
		  $temp2=substr($seq2,$sum2,$arr[$i-1]);  $new2="$new2$temp2";   $sum2=$sum2+$arr[$i-1];  
	   }
	 
	 $i++;
   }
   print DDD $line,"\t",$new1,"\t",$new2,"\n";
}
close CCC;  close DDD;

############
$len1=$small-1;
$addgene1=substr($ss1,0,$len1);  $addgene11="-"x$len1;
$len2=$genelen-$big;
$addgene2=substr($ss1,$big,$len2);  $addgene22="-"x$len2;
#print $big,"\t",$len2,"\t",$genelen,"\n";
############

$string1=""; $string2="";
$start1=0;  $start2=0;
$end1=0;  $end2=0;
open(CCC,"temp.out.sorted.result");
open(DDD,">temp.out.sorted.result2");
while($line=<CCC>)
{
  chomp $line;
  @bb=split/\s+/,$line;  
  if($bb[1]>$end1 && $bb[4]>$end2)
    {  $len1=abs($bb[1]-$end1)-1;
	   $len2=abs($bb[4]-$end2)-1;
	   $indel1=substr($ss1,$end1,$len1);  $add1="-"x$len1;  
	   $indel2=substr($ss2,$end2,$len2);  $add2="-"x$len2;  #print $bb[1],"\t",$end1,"\t",$len1,"\t",$bb[4],"\t",$end2,"\t",$len2,"\n";
	   if($end1==0 && $end2==0){$string1="$bb[7]";    $string2="$bb[8]";}
	   else{ $string1="$string1$indel1$add2$bb[7]";   $string2="$string2$add1$indel2$bb[8]";}
	}
    elsif($bb[1]<=$end1 && $bb[4]>$end2)
	{  $len1=abs($bb[1]-$end1)+1;
	   $len2=abs($bb[4]-$end2)-1;
	   $indel2=substr($ss2,$end2,$len2);  $add2="-"x$len2;
	   $string1="$string1$add2";  $string2="$string2$indel2";
	   @temparr=split//,$bb[7]; $len_temparr=@temparr; $ii=0;  $cc="-";
	   while($ii<$len1)
	   { $string1="$string1$cc";  $ii++;}
	   while($ii<$len_temparr)
	   { $string1="$string1$temparr[$ii]";  $ii++;}
	   $string2="$string2$bb[8]";
	}
    elsif($bb[1]>$end1 && $bb[4]<=$end2)
	{  $len1=abs($bb[1]-$end1)-1;
	   $len2=abs($bb[4]-$end2)+1;
	   $indel1=substr($ss1,$end1,$len1);  $add1="-"x$len1;
	   $string1="$string1$indel1";
	   $string1="$string1$bb[7]";
	   @temparr=split//,$bb[8]; $len_temparr=@temparr; $ii=0;
	   $string2="$string2$add1"; $cc="-";
	   while($ii<$len2)
	   { $string2="$string2$cc";  $ii++;}
	   while($ii<$len_temparr)
	   { $string2="$string2$temparr[$ii]";  $ii++;}
	}
	$start1=$bb[1];  $start2=$bb[4];
    $end1=$bb[2];    $end2=$bb[5];
}
print DDD  $addgene1,$string1,$addgene2,"\n",$addgene11,$string2,$addgene22,"\n";
close CCC;  close DDD;

open(DDD,"temp.out.sorted.result2");
$line1=<DDD>; $line2=<DDD>;
close DDD;
chomp $line1; chomp $line2;
$match="";
@a1=split//,$line1;
@a2=split//,$line2;
$len=@a1;
$i=0;
while($i<$len)
{
 if($a1[$i] eq $a2[$i])
 {
  $match="$match*";
 }
 else{$match="$match ";}
 $i++;
}
$temp=" "x30;

open(EEE,">bwasw.out.line");
print EEE  ">$scaffoldname\t$genename\t$direction\t$chr\t$small2\t$big2\n";
print EEE  "$genename$line1\n";
print EEE  "$scaffoldname$line2\n";
print EEE  "$temp$match\n";
close EEE;

}
#####################################
#####################################
#####################################
