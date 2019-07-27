#!/usr/bin/perl -w
#2017/4/25  yszhou

$bwain=$ARGV[0];  #"clustalw.out.checked.filtered";
$gene_structure=$ARGV[1];     #"gene.stucture.new";
$bwaout="$bwain.structure";   #"clustalw.out.checked.filtered.structure";

open(AA,$bwain);
while($s1=<AA>)
{ $s2=<AA>; $s3=<AA>; $s4=<AA>;
  $ll=$s1; chomp $ll; $ll=~s/>//;  
  @bb=split/\s+/,$ll;
  $protein=$bb[1]; $scaffold=$bb[0];
  $loc{$protein}=$scaffold;
  $ll2=$s2; $ll3=$s3; $ll4=$s4;
  chomp $ll2; chomp $ll3; chomp $ll4;
  @bb2=split//,$ll2; @bb3=split//,$ll3; @bb4=split//,$ll4;
  $i=0; $len=@bb2;
  while($bb2[$i] ne " "){$i++;}
  while($bb2[$i] eq " "){$i++;}
  $j=0;
  while($i<$len)
  {
  if($bb2[$i] ne "-"){
  $j++;
  $h2{"$protein,$j"}=$bb2[$i];
  $h3{"$protein,$j"}=$bb3[$i];
  $h4{"$protein,$j"}=$bb4[$i];
  }
  else{ 
  $temp2=$h2{"$protein,$j"};  $h2{"$protein,$j"}="$temp2$bb2[$i]";
  $temp3=$h3{"$protein,$j"};  $h3{"$protein,$j"}="$temp3$bb3[$i]";
  $temp4=$h4{"$protein,$j"};  $h4{"$protein,$j"}="$temp4$bb4[$i]";
  }
  $i++;
  }
}
close AA;

open(BB,"$gene_structure");
open(CC,">temp.stucture.utl-exon-intron");
while($line=<BB>)
{ chomp $line;
  if($line=~/^>/){$ll=$line; $ll=~s/>//; @bb=split/\s+/,$ll; $protein=$bb[0]; 
    if(exists($loc{$protein})){print CC $line,"\t",$loc{$protein},"\n";}
	}
  elsif(exists($loc{$protein})){
   @bb=split/\s+/,$line;
   $ss="$bb[0]  $bb[1]  $bb[2]";
   $len=30-length($ss); $temp=" "x$len;
   $i=$bb[1];
   $s2="$ss$temp"; $s3=" "x30;; $s4=" "x30;  
   while($i<=$bb[2])
   {  $t2=$h2{"$protein,$i"}; $t3=$h3{"$protein,$i"}; $t4=$h4{"$protein,$i"};
      $s2="$s2$t2";  $s3="$s3$t3";   $s4="$s4$t4";  
       $i++;
   }
     print CC $s2,"\n",$s3,"\n",$s4,"\n";
  }
}
close BB; close CC;

################################################################################################ Link the sequence
open(AA2,"temp.stucture.utl-exon-intron");
open(BB2,">temp.stucture.exon-intron.new");
$s=0;
while($line=<AA2>)
{  
  if($line=~/^>/)
    {if($s>0)
    {print BB2 $head,$s1,"\n".$s2,"\n",$s3,"\n",$s4,"\n",$s4,"\n";}
	$head=$line;
	$s1=""; $s2=""; $s3=""; $s4="";
    }
###########  
  if($line=~/^EXON/)
      {
	                chomp $line;  @b1=split//,$line;   
	  $line2=<AA2>;  chomp $line2; @b2=split//,$line2;  
	  $line3=<AA2>;  chomp $line3; @b3=split//,$line3;  
	  $len=@b1;
	  $i=30;  
	  while($i<$len)
	  {
	  $s1="$s1$b1[$i]";
	  $s2="$s2$b2[$i]";
	  $s3="$s3$b3[$i]";
	  $st="E";
	  $s4="$s4$st";
	  $i++;
	  }
	  }
##########
	  if($line=~/^INTRON/)
      {
	                chomp $line;  @b1=split//,$line;   
	  $line2=<AA2>;  chomp $line2; @b2=split//,$line2;  
	  $line3=<AA2>;  chomp $line3; @b3=split//,$line3;  
	  $len=@b1;
	  $i=30;  print $len,"\n";
	  while($i<$len)
	  {
	  $s1="$s1$b1[$i]";
	  $s2="$s2$b2[$i]";
	  $s3="$s3$b3[$i]";
	  $st="I";
	  $s4="$s4$st";
	  $i++;
	  }
	  }
###############	  
	 if($line=~/^UTL/)
      {
	                chomp $line;  @b1=split//,$line;   
	  $line2=<AA2>;  chomp $line2; @b2=split//,$line2;  
	  $line3=<AA2>;  chomp $line3; @b3=split//,$line3;  
	  $len=@b1;
	  $i=30;
	  while($i<$len)
	  {
	  $s1="$s1$b1[$i]";
	  $s2="$s2$b2[$i]";
	  $s3="$s3$b3[$i]";
	  $st="U";
	  $s4="$s4$st";
	  $i++;
	  }
	}
	 $s=1;
}
print BB2 $head,$s1,"\n".$s2,"\n",$s3,"\n",$s4,"\n",$s4,"\n";
close AA2; close BB2;
################################################################################################ Find ATG
open(AA3,"temp.stucture.exon-intron.new");
open(BB3,">temp.stucture.exon-intron.new.ATG");
while($s1=<AA3>)
{  
   $s2=<AA3>;  $s3=<AA3>; $s4=<AA3>;  $s5=<AA3>;  $s6=<AA3>;  
   
      $ll2=$s2; chomp $ll2;  @b2=split//,$ll2; $len=length($ll2);
      $ll3=$s3; chomp $ll3;  @b3=split//,$ll3;
      $ll5=$s5; chomp $ll5;  @b5=split//,$ll5;
      $ll6=$s6; chomp $ll6;  @b6=split//,$ll6;
      $i=0; 
	  while($b6[$i] ne "E" && $b6[$i] ne "F" ){$i++;}
      $j=$len-1;
	  while($b6[$j] ne "E" && $b6[$j] ne "F" ){$j--;}
      $start1="$b2[$i]$b2[$i+1]$b2[$i+2]";  $end1="$b2[$j-2]$b2[$j-1]$b2[$j]";
      $start2="$b3[$i]$b3[$i+1]$b3[$i+2]";  $end2="$b3[$j-2]$b3[$j-1]$b3[$j]";
	 
	 if($start1 eq "ATG" && $start2 ne "ATG" && $b6[$i] eq "E")
	  {
		  ####################### extend forward 300bp
		  $i2=$i+5;  $f="";
		  while($i2>200 && $f ne "ATG")
		  {
		  while($b3[$i2] eq "-"){$i2--;}
		  $t1=$b3[$i2];  $i2--;
		  while($b3[$i2] eq "-"){$i2--;}
		  $t2=$b3[$i2];  $i2--;
		  while($b3[$i2] eq "-"){$i2--;}
		  $t3=$b3[$i2];  $i2--;
		  $f="$t3$t2$t1";
		  }
		  if($f eq "ATG")
		  {$s7=""; $d="U"; $d2="E";  $i3=0; 
		  while($i3<=$i2){$s7="$s7$d";  $i3++;}
		  while($i3<=$i) {$s7="$s7$d2"; $i3++;}
		  while($i3<$len){$s7="$s7$b6[$i3]";  $i3++; }
		  chomp $s7; $s7="$s7\n";
		  #print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7;
		  }   ######################### extend 300+100bp
		 else{  
		   $i4=$i-300; $ru=0; $pre=""; $pre2="";
		   while($f ne "ATG" && $ru<400)
		   {
		      while($b3[$i4] eq "-"){$i4++;} 
			  $f="$pre2$pre$b3[$i4]";
			  $pre2=$pre; $pre="$b3[$i4]"; 
		      $i4++;  $ru++;
		   }
		  if($f eq "ATG" && $b6[$i4] ne "I")
		  {$s7=""; $d="U"; $d2="E";  $i5=0; 
		  while($i5<=($i4-4)){$s7="$s7$d";  $i5++;}
		  while($i5<=$i) {$s7="$s7$d2"; $i5++;}
		  while($i5<=$len){$s7="$s7$b6[$i5]";  $i5++; }
		  chomp $s7; $s7="$s7\n";
		  #print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7;
		  }################## can't find
		  else{$s7=$s6;}
	     }
	   }
	   else{$s7=$s6;}
	print BB3 $s1,$s2,$s3,$s4,$s5,$s6,$s7;
}
close AA3; close BB3;
################################################################################################ Find Stop Codon
open(AA4,"temp.stucture.exon-intron.new.ATG");
open(BB4,">temp.stucture.exon-intron.new.ATG.TAA"); 
while($s1=<AA4>)
{  
   $s2=<AA4>;  $s3=<AA4>;  $s4=<AA4>;  $s5=<AA4>;  $s6=<AA4>;  $s7=<AA4>;
   
      $ll2=$s2; chomp $ll2;  @b2=split//,$ll2; $len=length($ll2);
	  $ll3=$s3; chomp $ll3;  @b3=split//,$ll3;
	  $ll5=$s5; chomp $ll5;  @b5=split//,$ll5;
	  $ll6=$s6; chomp $ll6;  @b6=split//,$ll6;
	  $ll7=$s7; chomp $ll7;  @b7=split//,$ll7;
      $i=0; 
	  $exon=$ll7; $exon=~s/E//g; $len_exon=length($ll7)-length($exon);
	  if($len_exon>20){
	  while($b7[$i] ne "E"){$i++;}
      $j=$len-1;
	  while($b7[$j] ne "E"){$j--;}
      $start1="$b2[$i]$b2[$i+1]$b2[$i+2]";  $end1="$b2[$j-2]$b2[$j-1]$b2[$j]";
      $start2="$b3[$i]$b3[$i+1]$b3[$i+2]";  $end2="$b3[$j-2]$b3[$j-1]$b3[$j]";
      $ii=$i; $f=0; 
	  ###########
	  while($f==0 && $ii<$j)
	  { $pp=$ii;
	    while($b3[$ii] eq "-" || $b7[$ii] eq "I" || $b7[$ii] eq "F"){$ii++;}
		$a=$b3[$ii];  $ii++;
	    while($b3[$ii] eq "-" || $b7[$ii] eq "I" || $b7[$ii] eq "F"){$ii++;}
		$b=$b3[$ii];  $ii++;
		while($b3[$ii] eq "-" || $b7[$ii] eq "I" || $b7[$ii] eq "F"){$ii++;}
		$c=$b3[$ii];  
		$ii++;
		$temp="$a$b$c";  
		if($temp eq "TAA" || $temp eq "TAG" || $temp eq "TGA"){$f=1; }
	  }
       $ii=$pp;
      while($f==0 && $ii<($j+300))
	  { $pp=$ii;
	    while($b3[$ii] eq "-"){$ii++;}
		$a=$b3[$ii];  $ii++;
	    while($b3[$ii] eq "-"){$ii++;}
		$b=$b3[$ii];  $ii++;
		while($b3[$ii] eq "-"){$ii++;}
		$c=$b3[$ii];  
		$ii++;
		$temp="$a$b$c";  
		if($temp eq "TAA" || $temp eq "TAG" || $temp eq "TGA"){$f=1; }
	  } 
	  $s8="";
	  $ii2=0; $d="E"; $tt=0;
	  while($ii2<=$j) { $s8="$s8$b7[$ii2]"; $ii2++;}
	  while($ii2<$ii &&  $ii2<$len){ $s8="$s8$d";    $tt=1;    $ii2++;}
	  while($ii2<$len)
	  {$s8="$s8$b7[$ii2]"; $ii2++;}
	  chomp $s8; $s8="$s8\n";
	  print BB4 $s1,$s2,$s3,$s4,$s5,$s8;
	  }
}
close AA4; close BB4;

#######################################################################################  bwa2gene_structure-protein
$h{"TTT"}="F";  $h{"TTC"}="F";  $h{"TTA"}="L";  $h{"TTG"}="L";
$h{"TCT"}="S";  $h{"TCC"}="S";  $h{"TCA"}="S";  $h{"TCG"}="S";
$h{"TAT"}="Y";  $h{"TAC"}="Y";  $h{"TAA"}="*";  $h{"TAG"}="*";
$h{"TGT"}="C";  $h{"TGC"}="C";  $h{"TGA"}="*";  $h{"TGG"}="W";

$h{"CTT"}="L";  $h{"CTC"}="L";  $h{"CTA"}="L";  $h{"CTG"}="L";
$h{"CCT"}="P";  $h{"CCC"}="P";  $h{"CCA"}="P";  $h{"CCG"}="P";
$h{"CAT"}="H";  $h{"CAC"}="H";  $h{"CAA"}="Q";  $h{"CAG"}="Q";
$h{"CGT"}="R";  $h{"CGC"}="R";  $h{"CGA"}="R";  $h{"CGG"}="R";

$h{"ATT"}="I";  $h{"ATC"}="I";  $h{"ATA"}="I";  $h{"ATG"}="M";
$h{"ACT"}="T";  $h{"ACC"}="T";  $h{"ACA"}="T";  $h{"ACG"}="T";
$h{"AAT"}="N";  $h{"AAC"}="N";  $h{"AAA"}="K";  $h{"AAG"}="K";
$h{"AGT"}="S";  $h{"AGC"}="S";  $h{"AGA"}="R";  $h{"AGG"}="R";

$h{"GTT"}="V";  $h{"GTC"}="V";  $h{"GTA"}="V";  $h{"GTG"}="V";
$h{"GCT"}="A";  $h{"GCC"}="A";  $h{"GCA"}="A";  $h{"GCG"}="A";
$h{"GAT"}="D";  $h{"GAC"}="D";  $h{"GAA"}="E";  $h{"GAG"}="E";
$h{"GGT"}="G";  $h{"GGC"}="G";  $h{"GGA"}="G";  $h{"GGG"}="G";
#########################################

open(AA5,"temp.stucture.exon-intron.new.ATG.TAA");
open(BB5,">temp.stucture.exon-intron.new.ATG.TAA.protein");
while($s1=<AA5>)
{  $s2=<AA5>;
   $s3=<AA5>;
   $s4=<AA5>;
   $s5=<AA5>;
   $s6=<AA5>;
   ####
   $ll2=$s2; chomp $ll2;  @b2=split//,$ll2;  $len=@b2; $pre2=""; $pre1="";
   $ll5=$s5; chomp $ll5;  @b5=split//,$ll5;
   $ll6=$s6; chomp $ll6;  @b6=split//,$ll6;
   $i=0;  $t=0; $seq2=""; $stop=0;
   while($i<$len)
   {  if($b5[$i] eq "E" || $b5[$i] eq "F" )
      {
      if($b2[$i] ne "-" && $stop==0 )
      {$t++;
	  if($t<3){$pre2=$pre1; $pre1=$b2[$i]; $seq2="$seq2\_";}
        else{$temp="$pre2$pre1$b2[$i]"; 
		$seq2="$seq2$h{$temp}";
		if($h{$temp} eq "*"){$stop=1;}
		$pre2=$pre1; $pre1=$b2[$i];
	        $t=0;
			}
	   }
	  else{$seq2="$seq2 "}
	  }
	  else{$seq2="$seq2 ";}
	  $i++;
   }
   #####
   $ll3=$s3; chomp $ll3;  @b3=split//,$ll3;  $len=@b3; $pre2=""; $pre1="";
   $i=0;  $t=0; $seq3=""; $stop=0;
   while($i<$len)
   {  if($b6[$i] eq "E")
      {
      if($b3[$i] ne "-" && $stop==0)
      {$t++;
	  if($t<3){$pre2=$pre1; $pre1=$b3[$i]; $seq3="$seq3\_";}
        else{$temp="$pre2$pre1$b3[$i]";
		    $seq3="$seq3$h{$temp}";
			if($h{$temp} eq "*"){$stop=1;}
			$pre2=$pre1; $pre1=$b3[$i];
	        $t=0;
			}
	  }
	  else{$seq3="$seq3 "}
	  }
	  else{$seq3="$seq3 "}
	  $i++;
   }
     print BB5 $s1,$s2,$s3,$s4,$s5,$s6,$seq2,"\n",$seq3,"\n";
 }
 close AA5; close BB5;

##################################
`mv  temp.stucture.exon-intron.new.ATG.TAA.protein $bwaout`;
`rm  temp.stucture.exon-intron*`;  
