#!/usr/bin/perl
my $runs=0;
my $measures=0;

my @times;
my @vets;

my $datfile = $ARGV[0];

open(DAT,"<$datfile") || die " Could not open data file $datfile!\n";
my $line;
$line=<DAT>;
while( $line=<DAT>) {
    chomp($line);              # remove the newline from $line.

    my @pieces = split(/,/,$line);


    if ($pieces[0] eq "time") {
	$runs++;
	$measures=0;

    } else {
	if ($runs == 1) {
	    $times[$measures] = $pieces[0];
	}
	if ($pieces[1]==0) {
	    $vets[$measures] .= ",0";
	} else {
	    $vets[$measures] .= ",".log($pieces[1])/log(10);
	}
	$measures++;
    }
}
$runs++;
close(DAT);

printf ("time");
for (my $j=0; $j < $runs;$j++){
	printf (",vet%d",$j+1);
}
printf ("\n");
printf STDERR ("%d runs, %d measures\n",$runs,$measures);

for (my $j=0; $j < $measures;$j++){
    if ($times[$j] >= 150) {
	printf ("%g%s\n",$times[$j]-150.,$vets[$j]);
    }
}
