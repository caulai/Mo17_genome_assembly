#!/usr/bin/perl -w
#2017/4/20  yszhou

open(AA,"scaffold.fa");
$line1=<AA>;  
chomp $line1; $line1=~s/>//;  $chr=$line1; $len=length($line1);  $temp=" "x(30-$len);  $scaffoldname="$line1$temp";
$line2=<AA>;  
$ss2=$line2;  chomp $ss2;  
close AA;

####################
$find=0;  $good="";
$file="temp";
open(AA,"$file");
while($line=<AA>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if(@bb>8){
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
close AA;
`rm temp.sam2`;
if($good eq ""){`cp $file temp.sam2`;}
else{open(FA,">temp.sam2");  print FA $good,"\n"; close FA;}
#####################

$dir0=0; $dir16=0;
open(AA,"temp.sam2");
open(BB,">temp.out");
while($line=<AA>)
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
   print BB $bb[0],"\t",$start2,"\t",$end2,"\t",$bb[2],"\t",$start1,"\t",$end1,"\t",$bb[5],"\t",$bb[1],"\n";
  if($bb[1]==0){$dir0=$dir0+$match;}
    elsif($bb[1]==16){$dir16=$dir16+$match;}
  }
  }
}
close AA; close BB;
 if($dir0>$dir16){$direction="+"; $tt=0;}else{$direction="-"; $tt=16;}
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

$s=0;
open(AA,"temp.out.sorted");
open(BB,">temp.out.sorted.filterred");
while($line=<AA>)
{
  chomp $line;
  @bb=split/\s+/,$line;
  if($s>0){
    if($bb[4]<$start2)
     {if($prelen<($bb[2]-$bb[1])){$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}}
      elsif($bb[1]<$end1 && $bb[4]<$end2)
	   {if($prelen<($bb[2]-$bb[1])){$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}}
	  else{ print BB $preline,"\n"; $preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}
   }
   else{$preline=$line; $prelen=$bb[2]-$bb[1];  $start1=$bb[1];  $start2=$bb[4];  $end1=$bb[2];    $end2=$bb[5];}
   $s=1;
}
print BB $preline,"\n";
close AA; close BB;
###################################################################################################################################
if($direction eq "+"){$file="gene1.fa";}else{$file="gene2.fa";}
open(AA,$file);
$line1=<AA>;  chomp $line1; $line1=~s/>//;  $line1=~s/\r/ /; $len=length($line1);  $temp=" "x(30-$len);  $genename="$line1$temp";
$line2=<AA>;  $ss1=$line2;  chomp $ss1;    $genelen=length($ss1); 
close AA;

###################################################################################################################################

$small=10000000000;  $big=0;
$small2=10000000000;  $big2=0;
open(CC,"temp.out.sorted.filterred");
open(DD,">temp.out.sorted.result");
while($line=<CC>)
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
   print DD $line,"\t",$new1,"\t",$new2,"\n";
}
close CC;  close DD;
######
$len1=$small-1;
$addgene1=substr($ss1,0,$len1);  $addgene11="-"x$len1;
$len2=$genelen-$big;
$addgene2=substr($ss1,$big,$len2);  $addgene22="-"x$len2;
#print $big,"\t",$len2,"\t",$genelen,"\n";


#############
$sting1=""; $string2="";
$start1=0;  $start2=0;
$end1=0;  $end2=0;
open(CC,"temp.out.sorted.result");
open(DD,">temp.out.sorted.result2");
while($line=<CC>)
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
print DD  $addgene1,$string1,$addgene2,"\n",$addgene11,$string2,$addgene22,"\n";
close CC;  close DD;

open(DD,"temp.out.sorted.result2");
$line1=<DD>; $line2=<DD>;
close DD;
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
open(EE,">bwasw.out.line");
print EE  ">$scaffoldname\t$genename\t$direction\t$chr\t$small2\t$big2\n";
print EE  "$genename$line1\n";
print EE  "$scaffoldname$line2\n";
print EE  "$temp$match\n";
close EE;
