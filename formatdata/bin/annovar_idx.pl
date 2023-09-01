$N=0;
$bin=$ARGV[1];
my %idxMAX=();
my %idxMIN=();
open I,"<$ARGV[0]";
$KEY='XXX';
$idxMAX{'XXX'}=0;
while(<I>){
#	chomp;
	($a,$b)=(split/\t/)[0,1];
	$c=$bin*int($b/$bin);
	$length=length($_);
	$N+=$length;
	if($KEY eq "$a\t$c"){
		$idxMAX{$KEY}=$N;
	}else{
		$idxMIN{"$a\t$c"}=$idxMAX{$KEY};
		$KEY="$a\t$c";
		$idxMAX{$KEY}=$N;
	}
}
close I;

$size=`ls -l $ARGV[0]|awk \'{print \$5}\'`;
chomp($size);

print "#BIN\t$bin\t$size\n";
foreach(sort keys %idxMIN){
	if($_ eq 'XXX'){next;}
	if($_=~m/\#/){next;}
	print "$_\t$idxMIN{$_}\t$idxMAX{$_}\n";
}
