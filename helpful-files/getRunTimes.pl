use strict;
use warnings;
use Time::Piece;
#Nov-14 18:45:29.567 [Task submitter] INFO  nextflow.Session - [dc/a31060] Submitted process > impute (5)
#Nov-15 05:13:06.805 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 10; name: impute (5); status: COMPLETED; exit: 0; error: -; workDir: /pub40/harry/tgen/mgen/phase2/work/dc/a310608914b413aede3f6a312dc1a0]

my (%starts, %ends);
my $in = ".nextflow.log";
my $log = $ARGV[0];
if($log){
	$in = $log;
}
my @jobs = `grep "Submitted process" $in`;
foreach my $job (@jobs){
	my @parse = split(/\s+/, $job);
	if($parse[12]){
		$starts{"$parse[11] $parse[12]"}  = $parse[0] . " " . substr($parse[1], 0, -4);
	}else{
		$starts{$parse[11]}  = $parse[0] . " " . substr($parse[1], 0, -4);
	}
}

@jobs = `grep COMPLETED $in`;
foreach my $job (@jobs){
	my @parse = split(/\s+/, $job);
	if($parse[14] eq "status:"){
		$ends{$parse[13]}  = $parse[0] . " " . substr($parse[1], 0, -4);
	}else{
		$ends{"$parse[13] $parse[14]"}  = $parse[0] . " " . substr($parse[1], 0, -4);
	}
}
my $format = '%b-%d %H:%M:%S';
my %mean;
my %sd;
print"Process 	Start Time	End Time	Duration\n";
foreach my $k (sort keys  %starts){
	if($ends{"$k;"}){
		my $interval = Time::Piece->strptime( $ends{"$k;"}, $format) 
		- Time::Piece->strptime( $starts{$k}, $format);

		my ($proc, $no) = split(/\s+/,$k);
		$mean{$proc} += $interval;
		push(@{$sd{$proc}},$interval);
		print "$k\t$starts{$k}\t" . $ends{"$k;"} ."\t" . convert_seconds_to_hhmmss($interval). "\n";
	}
}
my $total_time;
foreach my $k (keys %mean){
	if($k ne "0"){
	my $me = 0;
	if(($#{$sd{$k}} + 1) . 0){
		$total_time += $mean{$k};
		$me = $mean{$k}/($#{$sd{$k}} + 1);
	}
		print "Mean time for process $k = " . convert_seconds_to_hhmmss($me) . "\n";
		print "Total time for process $k = " . convert_seconds_to_hhmmss($mean{$k}) . "\n";
	}
}
print "Total time for all processes = " . convert_seconds_to_hhmmss($total_time) . "\n";
##########################################
 sub convert_seconds_to_hhmmss {
  my $secs = shift;
  print "$secs";
  my $hourz=int($secs/3600);

  my $leftover=$secs % 3600;

  my $minz=int($leftover/60);

  my $secz=int($leftover % 60);

  

  return sprintf ("%02d:%02d:%02d", $hourz,$minz,$secz)

 

 }
