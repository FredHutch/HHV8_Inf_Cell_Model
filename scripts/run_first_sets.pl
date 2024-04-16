#!/usr/bin/perl
if ($#ARGV < 3) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 4 input arguments ($args given)\n";
}
my $infile;
my $outdir;
my $count;
my $model;

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
if ( $ARGV[3] =~ /^[12345678]$/) {
    $model=$ARGV[3]
} else {
    die "4th arg is the number of the model (1-7)\n";
    exit (1);
}
#an 0.5-2.5
my $low_an=0.5;
my $high_an=2.5;
my $range_an= $high_an-$low_an;

my $low_alpha=0.1;
my $high_alpha=0.5;
my $range_alpha= $high_alpha-$low_alpha;

my $low_log_r=1;
my $high_log_r=3;
my $range_log_r= $high_log_r-$low_log_r;

my $low_log_fpos=-2;
my $high_log_fpos=2;
my $range_log_fpos= $high_log_fpos-$low_log_fpos;

my $low_log_hill=-4;
my $high_log_hill=-1;
my $range_log_hill= $high_log_hill-$low_log_hill;

#Range for betaI (equiv to R0; NOT viral beta)
my $low_beta=1;
my $high_beta=20;
my $range_beta= $high_beta-$low_beta;

my $low_log_p=3.5;
my $high_log_p=7;
my $range_log_p= $high_log_p-$low_log_p;

my $low_log_r=1;
my $high_log_r=3;
my $range_log_r= $high_log_r-$low_log_r;

my $low_log_latent_inf=-2;
my $high_log_latent_inf=1;
my $range_log_latent_inf= $high_log_latent_inf-$low_log_latent_inf;

#my $low_expand_days=7;
#my $high_expand_days=20;
#my $range_expand_days= $high_expand_days-$low_expand_days;

# Models (1-5)
#1.Allow latency release rate & immune cell killing rate to vary
#2.Allow beta, latency release rate & immune cell killing rate to vary
#3.Allow latency release rate & infected cell death rate to vary
#4.Allow beta, latency release rate & infected cell death rate to vary
#5.Allow beta, immune cell killing rate, pi, r, latency release rate & infected cell death rate to vary

while ( $count > 0 ){
    print "$count runs remaining\n";

    my $beta= $range_beta*rand() + $low_beta;
    my $fpos= 10**($range_log_fpos*rand() + $low_log_fpos);
    my $hill= 10**($range_log_hill*rand() + $low_log_hill);
    my $log_p= $range_log_p*rand() + $low_log_p;
    my $r= 10**($range_log_r*rand() + $low_log_r);
    my $latent_inf= 10**($range_log_latent_inf*rand() + $low_log_latent_inf);
    my $alpha= $range_alpha*rand() + $low_alpha;
    my $an= $range_an*rand() + $low_an;
    #my $expand_days= $range_expand_days*rand() + $low_expand_days;

    open (TSTFILE, ">>$infile");
    printf TSTFILE ("PDF_on 0\n");
    if ( $model == 2 || $model == 4 || $model == 5 || $model == 8) {
	printf TSTFILE ("beta_init %e\n",$beta);
    }
    if ( $model == 1 || $model == 2 || $model >= 5) {
	printf TSTFILE ("fpos %g\n",$fpos);
    }
    if ( $model == 3 || $model == 4 || $model == 5 || $model == 8 ) {
	printf TSTFILE ("an %g\n",$an);
    }
    if ( $model < 6 ) {
	printf TSTFILE ("r_init %g\n",$r);
	printf TSTFILE ("alpha_init %g\n",$alpha);
    }
    if ( $model == 7 ) {
	printf TSTFILE ("density_killing 1\n");
	printf TSTFILE ("hill %g\n",$hill);
    }
    if ( $model == 8 ) {
	printf TSTFILE ("kappa_init 0\n");
    }
    printf TSTFILE ("log_p_init %g\n",$log_p);
    printf TSTFILE ("latent_inf_init %g\n",$latent_inf);
    #printf TSTFILE ("exp_days_init %g\n",$expand_days);
    close(TSTFILE);

    system("../hhv8_sim -b -r -f $infile -c $critfile");
    $count--;
}
