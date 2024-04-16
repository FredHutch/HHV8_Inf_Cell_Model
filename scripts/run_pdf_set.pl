#!/usr/bin/perl
if ($#ARGV < 2) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 3 input arguments ($args given)\n";
}
my $infile;
my $outdir;

if ( $ARGV[0] =~ /^[+-]?\d+$/) {
    $count=$ARGV[0]
} else {
    die "1st arg is the number of runs\n";
    exit (1);
}
if ( -f $ARGV[1]) {
    $infile=$ARGV[1]
} else {
    die "2nd arg is input file\n";
    exit (1);
}
if ( -f $ARGV[2]) {
    $critfile=$ARGV[2]
} else {
    die "3rd arg is the criteria file\n";
    exit (1);
}

while ( $count > 0 ){
    print "$count runs remaining\n";

    system("../hhv8_sim -b -r -f $infile -c $critfile -w 1");
    $count--;
}
