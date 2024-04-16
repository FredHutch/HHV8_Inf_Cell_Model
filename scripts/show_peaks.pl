#!/usr/bin/perl
use lib(".");
use dates;

my $num_ptids=0;
my @ptids;


my $ptid_file = $ARGV[0];
#ptid,arm,hiv,ks,spsite,collection_date,titer_ml,log10titer_ml,hhv8pos,enrolldt,hiv_vl_detect_enroll,DHC_,BegDate,EndDate,CD4,CD4_Date,RNA,RNA_Date,RNA_lcml,hiv_vl_detect

my $line;
my @swabs;
my $num_swabs=0;
my $num_ptids=0;
my $ptid_p=0;
my $ptid_first_date=0;

$swabs[$num_swabs]=0;

open(PTID,"<$ptid_file") || die " Could not open PTID file $ptid_file!\n";
my $line=<PTID>;
my $max_swab=0;
while( $line =<PTID>) {
    chomp($line);              # remove the newline from $line.

    my @pieces = split(/\,/,$line);
    my $ptid = int($pieces[0]);
    my $hhv8;

    if ($ptid != $ptid_p) {
	$num_ptids++;
	$ptid_first_date=$pieces[5];
	if ($ptid_p != 0) {
	    my $crit_file="ptid_".$ptid_p.".crit";
	    if (open(CRIT,"<$crit_file") ){
		$line=<CRIT>;
		close(CRIT);

		chomp($line);              # remove the newline from $line.
		$line.=" ".$max_swab;
		open(CRIT,">$crit_file") || die " Could not open PTID crit file $crit_file!\n";
		printf CRIT ("%s\n",$line);
		close(CRIT);
		printf STDERR ("PTID %d: max swab level = %g\n",$ptid_p,$max_swab);
	    }
	}
	$max_swab=0;
    } else {
	if (sub_dates($pieces[5],$ptid_first_date) > 28) {
	    next;
	}
    }
    if ($pieces[8] != "") {
	$hhv8 = int($pieces[8]);
    } else {
	$hhv8 = 0;
    }
    if ($pieces[7] != "") {
	$hhv8_level = $pieces[7];
	if ($hhv8_level > $max_swab) {
	    $max_swab=$hhv8_level;
	}
    } else {
	$hhv8_level = 0;
    }
    $swabs[$num_swabs]=$hhv8_level;
    $num_swabs++;
    $ptid_p=$ptid;
}
close(PTID);

my $crit_file="ptid_".$ptid_p.".crit";
if (open(CRIT,"<$crit_file")){
    $line=<CRIT>;
    close(CRIT);

    chomp($line);              # remove the newline from $line.
    $line.=" ".$max_swab;
    open(CRIT,">$crit_file") || die " Could not open PTID crit file $crit_file!\n";
    printf CRIT ("%s\n",$line);
    close(CRIT);

    printf STDERR ("PTID %d: max swab level = %g\n",$ptid_p,$max_swab);
}

print STDERR "Read $num_ptids participants from $ptid_file ($num_swabs swabs)\n";
