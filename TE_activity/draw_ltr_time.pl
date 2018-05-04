use SVG;
use strict;
die "USAGE: perl B73.genome.len Mo17.genome.len B73.ltr.time Mo17.ltr.time\n" if(@ARGV!=4);
my %len_b73;
my %len_mo17;
my %count_b73;
my %count_mo17;

open A,$ARGV[0] or die $!;
while(<A>)
{
	chomp;
	my @a=split;
	$len_b73{$a[0]}=$a[1];
	for(my $i=0;$i<$a[1];$i+=1000000)
	{
		$count_b73{$a[0]}{$i}=0;
	}
}
close A;

open B,$ARGV[1] or die $!;
while(<B>)
{
	chomp;
	my @a=split;
	$len_mo17{$a[0]}=$a[1];
	for(my $i=0;$i<$a[1];$i+=1000000)
	{
		$count_mo17{$a[0]}{$i}=0;
	}
}
close B;

my $svg=SVG->new(width=>1850,height=>2000);
my $bp_per_pix=200000;
my $x_start=175;
my $y_start=0;
my %pos_b73;
my %pos_mo17;
for(my $i=1;$i<11;$i++)
{
	my $x1=$x_start;
	my $y1=($i-1)*180+100;
	my $chr="chr".$i;
	my $x2=$len_b73{$chr}/$bp_per_pix + $x_start;
	my $y2=($i-1)*180+100;
	$pos_b73{$chr}=$y1;
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke'=>'white','stroke-width'=>50,'fill'=>'white');
	my $chr_name1="B73_".$chr;
	my $text_x=20;
	my $text_y=$y1+6;
	$svg->text(x=>$text_x,y=>$text_y,'font-size'=>25,)->cdata("$chr_name1");
	my $x1=$x_start;
	my $y1=($i-1)*180+170;
	my $x2=$len_mo17{$chr}/$bp_per_pix + $x_start;
	my $y2=($i-1)*180+170;
	$pos_mo17{$chr}=$y1;
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'white','stroke-width'=>50,fill=>'white');
	my $chr_name2="Mo17_".$chr;
	my $text_x=15;
	my $text_y=$y1+6;
	$svg->text(x=>$text_x,y=>$text_y,'font-size'=>25,)->cdata("$chr_name2");
}

my $kedu_y=1900;
my $kedu=int($len_b73{chr1}/10000000)+1;
my $x1=$x_start;
my $y1=$kedu_y;
my $x2=$x_start+$kedu*10000000/$bp_per_pix;
my $y2=$kedu_y;
$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke','black','stroke-width',3,'fill','black',);
for(my $i=0;$i<=$kedu;$i++)
{
	if($i%5==0)
	{
		my $x1=$x_start+$i*10000000/$bp_per_pix;
		my $y1=$y_start+$kedu_y-20;
		my $x2=$x_start+$i*10000000/$bp_per_pix;
		my $y2=$y_start+$kedu_y;
		$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke','black','stroke-width',3,'fill','black',);
		my $text_x;
		my $text_y;
		my $cdata;
		if($i==0)
		{
			$cdata=0;
			$text_x=$x_start+$i*10000000/$bp_per_pix-4;
			$text_y=$y_start+$kedu_y+35;
		}
		else
		{
			$cdata=$i*10;$cdata=$cdata."Mb";
			$text_x=$x_start+$i*10000000/$bp_per_pix-25;
			$text_y=$y_start+$kedu_y+35;
		}
		$svg->text(x=>$text_x,y=>$text_y-2,'font-size'=>25,)->cdata("$cdata");
	}
	else
	{
		my $x1=$x_start+$i*10000000/$bp_per_pix;
		my $y1=$y_start+$kedu_y-10;
		my $x2=$x_start+$i*10000000/$bp_per_pix;
		my $y2=$y_start+$kedu_y;
		$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke','black','stroke-width',2,'fill','black',);
	}
}

my %num_b73;
open C,$ARGV[2] or die $!;
while(<C>)
{
	chomp;
	next if(/^>B73/);
	my @a=split(/\t/);
	$a[0]=~s/>//;
	my @b=split(/_/,$a[0]);
	my $chr_1="chr".$b[0];
	my $i=(int($b[-2]/1000000))*1000000;
	$count_b73{$chr_1}{$i}+=$a[-1];
	$num_b73{$chr_1}{$i}+=1;
}
close C;

open TEMP1,">$ARGV[2].count" or die $!;
foreach my $key1 (keys %count_b73)
{
	foreach my $key2 (keys %{$count_b73{$key1}})
	{
		my $end=$key2+1000000;
		my $value;
		if($num_b73{$key1}{$key2}==0){$value="un";}
		else{$value=$count_b73{$key1}{$key2}/$num_b73{$key1}{$key2};}
		print TEMP1 "$key1\t$key2\t$end\t$value\n";
	}
}
close TEMP1;

open C,"$ARGV[2].count" or die $!;
while(<C>)
{
	chomp;
	next if(/^>B73/);
	my @a=split(/\t/);
	my $x1=$a[1]/$bp_per_pix+$x_start;
	my $y1=$pos_b73{$a[0]};
	my $x2=$a[2]/$bp_per_pix+$x_start;
	my $y2=$pos_b73{$a[0]};
	if($a[-1]>1.5){$a[-1]=1.5;}
	my $color1=255-int($a[-1]*170);
	my $color2=int($a[-1]*170);
	my $fill="rgb($color1,0,$color2)";
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke',"$fill",'stroke-width',50,'fill',"$fill",);
}
close C;

my %num_mo17;
open D,$ARGV[3] or die $!;
while(<D>)
{
	chomp;
	next if(/^>scaffold/);
	my @a=split(/\t/);
	$a[0]=~s/>//;
	my @b=split(/_/,$a[0]);
	my $i=(int($b[-2]/1000000))*1000000;
	$count_mo17{$b[0]}{$i}+=$a[-1];
	$num_mo17{$b[0]}{$i}+=1;
}
close D;

open TEMP2,">$ARGV[3].count" or die $!;
foreach my $key1 (keys %count_mo17)
{
	foreach my $key2 (keys %{$count_mo17{$key1}})
	{
		my $end=$key2+1000000;
		my $value;
		if($num_mo17{$key1}{$key2}==0){$value="un";}
		else{$value=$count_mo17{$key1}{$key2}/$num_mo17{$key1}{$key2};}
		print TEMP2 "$key1\t$key2\t$end\t$value\n";
        }
}
close TEMP2;

open D,"$ARGV[3].count" or die $!;
while(<D>)
{
	chomp;
	my @a=split(/\t/);
	my $x1=$a[1]/$bp_per_pix+$x_start;
	my $y1=$pos_mo17{$a[0]};
	my $x2=$a[2]/$bp_per_pix+$x_start;
	my $y2=$pos_mo17{$a[0]};
	if($a[-1]>1.5){$a[-1]=1.5;}
	my $color1=255-int($a[-1]*170);
	my $color2=int($a[-1]*170);
	my $fill="rgb($color1,0,$color2)";
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke',"$fill",'stroke-width',50,'fill',"$fill",);
}
close D;

open E,"/NAS7/home/sunsl/mo17_analysis/b73_v4/centromere_maize_v4.txt" or die $!;
while(<E>)
{
	chomp;
	my @a=split(/\t/);
	my $x1=$a[1]/$bp_per_pix+$x_start;
	my $y1=$pos_b73{$a[0]}-35;
	my $x2=$a[2]/$bp_per_pix+$x_start;
	my $y2=$pos_b73{$a[0]}-35;
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke','black','stroke-width',20,'fill','black',);
}
close E;

for(my $i=0;$i<=255;$i+=5)
{
	my $color1=255-$i;
	my $color2=$i;
	my $fill="rgb($color1,0,$color2)";
	my $x1=1450+$i;
	my $y1=1700;
	my $x2=1450+$i+5;
	my $y2=1700;
	$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,'stroke',"$fill",'stroke-width',20,'fill',"$fill",);
	$svg->line(x1=>$x1,y1=>$y1+40,x2=>$x2,y2=>$y2+40,'stroke',"black",'stroke-width',5,'fill','black',);
	if($i%40==0)
	{
		$svg->line(x1=>$x1+2,y1=>$y1+40,x2=>$x1+2,y2=>$y2+50,'stroke',"black",'stroke-width',3,'fill','black',);
		my $cdata=$i/160;
		$svg->text(x=>$x1-10,y=>$y1+70,'font-size'=>20,)->cdata("$cdata");
	}
}

my $svg_name="B73_Mo17_LTR_date.svg";
open OUTPUT, ">$svg_name";
print OUTPUT $svg->xmlify;
close OUTPUT;

system("/NAS2/sunsl/pav/distributing_svg_4.74/svg2xxx_release/svg2xxx ".$svg_name." -t png");
