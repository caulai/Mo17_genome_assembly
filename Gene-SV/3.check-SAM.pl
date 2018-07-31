#!/usr/bin/perl -w
#2017/4/22  

$ref_genome=$ARGV[0];  # the second genome  (Mo17/B73)
$cds_seq=$ARGV[1];     # "gene-full-cds-double.new.fa"

################################################################################################################
open(FF,"$cds_seq");
while($line=<FF>)
{ chomp $line;
  if($line=~/^>/)
  { $ll=$line; $ll=~s/>//;  @bb=split/\s+/,$ll; $pre=$bb[0];  $dir=$bb[1];}
  else{
   if($dir==1){$ge{$pre}=$line;}
   else{$ge2{$pre}=$line;}
    }
}
close FF;
################################################################################################################
open(AA,"$ref_genome");
while($line1=<AA>)
{  $line2=<AA>;
   chomp $line1; chomp $line2; 
   $line1=~s/>//;  $line1=~s/chr0//;  $line1=~s/chr//; $line1=~s/Chr//;
   @bb=split/\s+/,$line1;
   $seq{$bb[0]}=$line2;  
   #print $bb[0],"\n";
}close AA;
################################################################################################################
open(AA,"bwasw.out");
open(BB,">bwasw.out.checked");
open(CC,">bad");
while($s1=<AA>)
{
  $s2=<AA>;  $s3=<AA>; $s4=<AA>;
  $ll1=$s1;  chomp $ll1;  @b1=split/\s+/,$ll1;  $start=$b1[4];  $end=$b1[5];  $mys=$b1[3];  $mys=~s/chr0//;  $mys=~s/chr//; $mys=~s/Chr//;
  $ll2=$s2;  chomp $ll2;  @b2=split/\s+/,$ll2;  $temp=$b2[1];   $temp=~s/-//g;
  $ll3=$s3;  chomp $ll3;  @b3=split/\s+/,$ll3;  $temp3=$b3[1];  $temp3=~s/-//g;
  $ll4=$s4;  chomp $ll4;  $temp2=$ll4;  $temp22=$ll4;  $temp22=~s/\*//g;
  ######
  $len1=$end-$start+1;  
  if($len1<1000000 && $len1>10)
  {  
  $st1=substr($seq{$mys},$start-1,$len1);  
  print $len1,"\t",length($temp3),"\t",length($st1),"\t",length($temp),"\t",length($ge{$b1[1]}),"\t",length($ge2{$b1[1]}),"\n"; 
  $st2=$st1;
  $st2=~tr/ATCGatcg/TAGCtagc/;
  @ar=split//,$st2;
  $lenar=@ar;
  $i=0;  $st2="";
  while($i<$lenar){$st2="$ar[$i]$st2"; $i++;}
  ######
  if($temp eq $ge{$b1[1]}  ||  $temp eq $ge2{$b1[1]})
   { if($temp3 eq $st1  || $temp3 eq $st2)
    { 
   if((length($temp2)-length($temp22))>100)
    { 
    ###############
    if($temp eq $ge{$b1[1]})
    {print BB $s1,$s2,$s3,$s4;}
    else
    {
     print BB $s1;
     @arr2=split//,$ll2;  	 @arr3=split//,$ll3;    @arr4=split//,$ll4; 
	   $len=@arr2;
	   $i=0; $st2=""; $st3=""; $st4="";
	   while($i<30)
	  {$st2="$st2$arr2[$i]"; $st3="$st3$arr3[$i]"; $st4="$st4$arr4[$i]"; $i++;}
	  $j=$len;
	  while($j>=$i)
	   { $t2=$arr2[$j-1];  $t3=$arr3[$j-1];  $t4=$arr4[$j-1];
	    $t2=~s/A/x/; $t2=~s/T/A/; $t2=~s/x/T/;   $t2=~s/G/y/; $t2=~s/C/G/; $t2=~s/y/C/; 
		  $t3=~s/A/x/; $t3=~s/T/A/; $t3=~s/x/T/;   $t3=~s/G/y/; $t3=~s/C/G/; $t3=~s/y/C/; 
	    $st2="$st2$t2"; $st3="$st3$t3"; $st4="$st4$t4"; 
     	$j--;
	   }
	   print BB $st2,"\n",$st3,"\n",$st4,"\n";
	   }  
     print "find the shit\n";
     ##############
     }
   }
   else{print CC $b1[1],"\t","scaffold","\n";}
  }
  else{print CC $b1[1],"\t","gene","\n";}
 }
}
close AA; close BB;  close CC;

########################################################################################################
open(AA,"bwasw.out.checked");
open(BB,">clustalw.out.match.dir.lenth");
while($s1=<AA>)
{  chomp $s1;  @b=split/\s+/,$s1;
   $s2=<AA>; chomp $s2;  @b2=split//,$s2; $len=@b2;
   $s3=<AA>; chomp $s3;  @b3=split//,$s2;
   $s4=<AA>; chomp $s4;  @b4=split//,$s2;
   $i=30; $j=0;
   while($i<$len)
   { if($b2[$i] ne "-"){$j++;}
     if($j==2000){$gstart=$i;}
	 $i++;
   }
   $i=$len; $j=0;
   while($i>30)
   { if($b2[$i] ne "-"){$j++;}
     if($j==2000){$gend=$i;}
     $i--;
   }
   $matchlen=$gend-$gstart;
   $matchcds=substr($s4,$gstart,$matchlen);
   $matchcds2=$matchcds;  $matchcds2=~s/\*//g;
   $matchgene=$s4; $matchgene=~s/\*//g;
   $lenth1=length($matchcds)-length($matchcds2);
   $lenth2=length($s4)-length($matchgene);
   print BB "$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$lenth1\t$lenth2\n";
}
close AA; close BB;
`sort -dk2,2 -k6,6nr  -k7,7nr  clustalw.out.match.dir.lenth > clustalw.out.match.dir.lenth.sorted`;
open(AA,"clustalw.out.match.dir.lenth.sorted");
open(BB,">clustalw.out.match.dir.lenth.sorted.filtered");
while($line=<AA>)
{ chomp $line; 
  @bb=split/\s+/,$line;
  if(not exists($hash{$bb[1]}))
  { $hash{$bb[1]}=1; 
    print BB $line,"\n"; $h{"$bb[0],$bb[1],$bb[2],$bb[3],$bb[4]"}=1;
  }
}
close AA;
##########################################
open(AA,"bwasw.out.checked");
open(BB,">clustalw.out.checked.filtered");
while($s1=<AA>)
{  $s2=<AA>;
   $s3=<AA>;
   $s4=<AA>;
   $ll=$s1;  chomp $ll; @bb=split/\s+/,$ll;
   if(exists($h{"$bb[0],$bb[1],$bb[2],$bb[3],$bb[4]"})){print BB $s1,$s2,$s3,$s4;} 

}
close AA; close BB;
`rm  clustalw.out.match.dir.lenth`;
`rm  clustalw.out.match.dir.lenth.sorted`;
`rm  clustalw.out.match.dir.lenth.sorted.filtered`;
##################################################################################################

