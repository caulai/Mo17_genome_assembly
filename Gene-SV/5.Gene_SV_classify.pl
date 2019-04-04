#!/usr/bin/perl -w
#2017/4/27  yszhou

$filein=$ARGV[0];  #clustalw.out.checked.filtered.structure

open(AA,"$filein");
open(BB,">$filein.large.variation.class3.CDS.incomplete");
open(CC,">$filein.clas12.CDS.complete");
while($s1=<AA>)
{  $s2=<AA>; $ll2=$s2;  #chomp $ll2; @b2=split//,$ll2;  
   $s3=<AA>; $ll3=$s3;  chomp $ll3; @b3=split//,$ll3;  
   $s4=<AA>; #$ll4=$s4;  chomp $ll4; @b4=split//,$ll4;  
   $s5=<AA>; $ll5=$s5;  chomp $ll5; @b5=split//,$ll5; 
   $s6=<AA>;
   $s7=<AA>;
   $s8=<AA>;
   $len=@b5;  $find=0;
   $i=0;  #$match=0;  $all=0;
   while($i<$len)
   {
     if($b5[$i-1] ne "E"  && $b5[$i] eq "E")
     { #$temp="$b3[$i-3]$b3[$i-2]$b3[$i-1]$b3[$i]$b3[$i+1]$b3[$i+2]";
       $temp="$b3[$i]$b3[$i+1]$b3[$i+2]";
       if($temp eq "---"){$find=1;}
     }
     if($b5[$i] eq "E"  && $b5[$i+1] ne "E")
     { #$temp="$b3[$i-2]$b3[$i-1]$b3[$i]$b3[$i+1]$b3[$i+2]$b3[$i+3]";
       $temp="$b3[$i-2]$b3[$i-1]$b3[$i]";
       if($temp eq "---"){$find=1;}
     }
       $i++;
   }
   if($find>0){print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
   else{print CC $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
}
close AA; close BB; close CC;
#########################################################################################################
open(AA,"$filein.clas12.CDS.complete");
open(F1,">$filein.class2.with_3n+-1_indel_in_CDS");
open(F2,">$filein.class2.Start_condon_mutation");
open(F3,">$filein.class2.Stop_codon_mutation");
open(F4,">$filein.class2.Splice_donor_mutation");
open(F5,">$filein.class2.Splice_acceptor_mutation");
open(F6,">$filein.class2.Premature_stop_codon");
open(F7,">$filein.class1.Without_amino_acid_substitutions");
open(F8,">$filein.class1.With_amino_acid_changes");
open(F9,">$filein.class2_Gene_with_Large_effect_mutations");
open(F10,">$filein.class1.No_DNA_variation_in_genetic_region");
open(F11,">$filein.class1.No_DNA_variation_in_CDS_intron_region");
open(F12,">$filein.class1.Structurally_conserved_genes");
###########################################################################################################
while($s1=<AA>)
{  $s2=<AA>;  $ll2=$s2;  chomp $ll2; @b2=split//,$ll2;   #seq1
   $s3=<AA>;  $ll3=$s3;  chomp $ll3; @b3=split//,$ll3;   #seq2
   $s4=<AA>;  $ll4=$s4;  chomp $ll4; @b4=split//,$ll3;
   $s5=<AA>;  $ll5=$s5;  chomp $ll5; @b5=split//,$ll5;   #structure
   $s6=<AA>;
   $s7=<AA>;  $ll7=$s7;  chomp $ll7; @b7=split//,$ll7;   #protein1
   $s8=<AA>;  $ll8=$s8;  chomp $ll8; @b8=split//,$ll8;   #protein2
   $f=0;
   $i=0; $len=@b7;  $st2="";   $cdsintron1="";    $cdsintron2="";  #$p1=""; $p2="";  # $st1="";
   $Findel=0;
   $FATG=0;
   $Fstop=0;
   $Fdonor=0;
   $Facceptor=0;
   $Fpremature=0;
   while($i<$len)
   {  #######################################################################################CDS-intron
      if($b5[$i] eq "E" || $b5[$i] eq "I"){$cdsintron1="$cdsintron1$b2[$i]";    $cdsintron2="$cdsintron2$b3[$i]";}
      ####################################################################################### indel
      if($b5[$i] eq "E" )
	  {  if($b3[$i] eq "-" && $b3[$i-1] ne "-"){$indel1=1;}
	     elsif($b3[$i] eq "-" && $b3[$i-1] eq "-"){$indel1++;}
		 elsif($b3[$i] ne "-" && $b3[$i-1] eq "-"){ $rate=$indel1%3; if($rate>0){$Findel=1;}}
	  
	     if($b4[$i] eq "-" && $b4[$i-1] ne "-"){$indel2=1;}
		 elsif($b4[$i] eq "-" && $b4[$i-1] eq "-"){$indel2++;}
	     elsif($b4[$i] ne "-" && $b4[$i-1] eq "-"){ $rate=$indel2%3; if($rate>0){$Findel=1;}}
	  }
	  ######################################################################################## ATG
	  if($b5[$i] eq "E" && $b5[$i-1] eq "U")
	  {  if($b2[$i] ne "-" && $b2[$i+1] ne "-" && $b2[$i+2] ne "-" && $b3[$i] ne "-" && $b3[$i+1] ne "-" && $b3[$i+2] ne "-")
	      {  $stp1="$b2[$i]$b2[$i+1]$b2[$i+2]";  $stp2="$b3[$i]$b3[$i+1]$b3[$i+2]"; 
	         if($stp1 ne $stp2){$FATG=1;}
	      }
	  }
	  ######################################################################################## StopCodon
	  if($b5[$i] eq "E" && $b5[$i+1] eq "U")
	  {  if($b2[$i-2] ne "-" && $b2[$i-1] ne "-" && $b2[$i] ne "-" && $b3[$i-2] ne "-" && $b3[$i-1] ne "-" && $b3[$i] ne "-")
	      { $stp1="$b2[$i-2]$b2[$i-1]$b2[$i]";  $stp2="$b3[$i-2]$b3[$i-1]$b3[$i]"; 
	        if($stp2 eq "TAA" || $stp2 eq "TGA" || $stp2 eq "TAG"){$Fstop=0;}
		    else{$Fstop=1;}
		    if($stp1 eq $stp2){$Fstop=0;}
		  }
	  }
	  ######################################################################################## donor & acceptor
	  if($b5[$i] eq "E" && $b5[$i+1] eq "I")
	  { if($b2[$i+1] ne "-" && $b2[$i+2] ne "-" && $b3[$i+1] ne "-" && $b3[$i+2] ne "-")
	    { $sp1="$b2[$i+1]$b2[$i+2]";  $sp2="$b3[$i+1]$b3[$i+2]";
	      if($sp1 ne $sp2){if($sp2 ne "GC" && $sp2 ne "GT"){$Fdonor=1;}}
	    }
	  }
	  if($b5[$i] eq "E" && $b5[$i-1] eq "I")
	  { if($b2[$i-2] ne "-" && $b2[$i-1] ne "-" && $b3[$i-2] ne "-" && $b3[$i-1] ne "-")
	    { $sp1="$b2[$i-2]$b2[$i-1]";  $sp2="$b3[$i-2]$b3[$i-1]";
	      if($sp1 ne $sp2){$Facceptor=1;}
	    }
	  }
	  ######################################################################################## premature
	  if($b5[$i] eq "E" && $b5[$i+1] eq "E" && $b8[$i] eq "*")
	  {   if($b2[$i] ne "-" && $b2[$i-1] ne "-" && $b2[$i-2] ne "-" && $b3[$i] ne "-" && $b3[$i-1] ne "-" && $b3[$i-2] ne "-")
	      {  $stp1="$b2[$i-2]$b2[$i-1]$b2[$i]";   $stp2="$b3[$i-2]$b3[$i-1]$b3[$i]"; 
		     if($stp1 ne $st2){$Fpremature=1;}
	      }
	  }
	  ########################################################################################
	$i++;
   } 
    if($Findel>0){print F1 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($FATG>0){print F2 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($Fstop>0){print F3 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($Fdonor>0){print F4 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($Facceptor>0){print F5 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($Fpremature>0){print F6 $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($Findel==0 && $FATG==0 && $Fstop==0 && $Fdonor==0 && $Facceptor==0 && $Fpremature==0)
	{ if($s7 eq $s8){print F7  $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	  else{print F8  $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
           print F12   $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;
	}
	else{print F9  $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($s2 eq $s3){print F10  $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
	if($cdsintron1 eq $cdsintron2){print F11  $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
} 
close AA; close F1; close F2; close F3; close F4; close F5; close F6; close F7; close F8; close F9; close F10; close F11;  close F12;
###############################################################################################################
open(AA,"$filein.class1.Without_amino_acid_substitutions");
open(BB,">$filein.class1.No_DNA_variation_in_CDS_region");
while($s1=<AA>)
{  $s2=<AA>;  $ll2=$s2;  chomp $ll2; @b2=split//,$ll2;   #seq1
   $s3=<AA>;  $ll3=$s3;  chomp $ll3; @b3=split//,$ll3;   #seq2
   $s4=<AA>;
   $s5=<AA>;  $ll5=$s5;  chomp $ll5; @b5=split//,$ll5;   #structure
   $s6=<AA>;
   $s7=<AA>;  $ll7=$s7;  chomp $ll7; @b7=split//,$ll7;   #protein1
   $s8=<AA>;  $ll8=$s8;  chomp $ll8; @b8=split//,$ll8;   #protein2
   $len=@b7;
   $i=0; $ss1=""; $ss2="";
   while($i<$len)
   { if($b5[$i] eq "E"){$ss1="$ss1$b2[$i]"; $ss2="$ss2$b3[$i]";} $i++;}
   if($ss1 eq $ss2){print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
   
}
close AA; close BB;
#################################################################################################################
open(AA,"$filein.large.variation.class3.CDS.incomplete");
open(BB,">$filein.large.variation.class3.CDS.incomplete.No_Exon_missing");
while($s1=<AA>)
{  $s2=<AA>;  $ll2=$s2;  chomp $ll2; @b2=split//,$ll2;   #seq1
   $s3=<AA>;  $ll3=$s3;  chomp $ll3; @b3=split//,$ll3;   #seq2
   $s4=<AA>;
   $s5=<AA>;  $ll5=$s5;  chomp $ll5; @b5=split//,$ll5;   #structure
   $s6=<AA>;
   $s7=<AA>;  $ll7=$s7;  chomp $ll7; @b7=split//,$ll7;   #protein1
   $s8=<AA>;  $ll8=$s8;  chomp $ll8; @b8=split//,$ll8;   #protein2
   $len=@b7;
   $i=0; $f=0; $find=0;
   while($i<$len)
   { if($b5[$i] eq "E" && $b5[$i-1] ne "E"){if($b3[$i] eq "-"){$f=1;}else{$f=0;} }
	 if($b5[$i] eq "E" && $b5[$i-1] eq "E")
	 {  if($b3[$i] ne "-"){$f=0;}}
	 if($b5[$i] eq "E" && $b5[$i+1] ne "E")
	 {  if($b3[$i] eq "-" && $f>0){$find=1;}}
	 $i++;
   }# print $find,"\n";
   if($find==0){print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
}
close AA; close BB;

###############################################################
open(AA,"$filein.class1.With_amino_acid_changes");
open(BB,">$filein.class1.With_3n_indel_in_CDS");
open(CC,">$filein.class1.With_missence_mutation_in_CDS");
while($s1=<AA>)
{  $s2=<AA>;  $ll2=$s2;  chomp $ll2; @b2=split//,$ll2;   #seq1
   $s3=<AA>;  $ll3=$s3;  chomp $ll3; @b3=split//,$ll3;   #seq2
   $s4=<AA>;
   $s5=<AA>;  $ll5=$s5;  chomp $ll5; @b5=split//,$ll5;   #structure
   $s6=<AA>;
   $s7=<AA>;  $ll7=$s7;  chomp $ll7; @b7=split//,$ll7;   #protein1
   $s8=<AA>;  $ll8=$s8;  chomp $ll8; @b8=split//,$ll8;   #protein2
   $len=@b7;
   $i=0; $ss1=""; $ss2="";
   $indel=0; $mis=0;
   while($i<$len)
   { if($b5[$i] eq "E"){
     if($b7[$i] eq "_" && $b8[$i] ne "_"){$indel=1;}
     if($b7[$i] ne "_" && $b8[$i] eq "_"){$indel=1;}
     if($b7[$i] ne "_" && $b8[$i] ne "_" && $b7[$i] ne $b8[$i]){$mis=1;}
     }
    $i++;
   }
   if($indel>0){print BB $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
   if($mis>0){print CC $s1,$s2,$s3,$s4,$s5,$s6,$s7,$s8;}
}

close AA; close BB;  close CC;
