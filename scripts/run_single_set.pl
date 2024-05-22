#!/usr/bin/perl
# Function to find random numbers that follow a normal distribution
# $ave is mean, $sigma is standard deviation,
sub randn {
  my ($ave, $sigma) = @_;
  my $ret=0;
  
  printf STDERR ("for mean=%g, stddev=%g..\n",$ave, $sigma);
  do {
      my ($r1, $r2) = (rand(), rand());
      while ($r1 == 0) {$r1 = rand();}
      $ret = ($sigma * sqrt(-2 * log($r1)) * sin(2 * 3.14159265359 * $r2)) + $ave;
  } while ($ret < 0);
  printf STDERR ("gaussian draw=%g\n",$ret);
  return $ret;
}

if ($#ARGV < 3) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 4 input arguments ($args given)\n";
}
my $infile;
my $outdir;
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

if ( $ARGV[3] =~ /^[123]/) {
    $model=$ARGV[3]
} else {
    die "4th arg is the number of the model (1-3)\n";
    exit (1);
}

my $mean_an;
my $stddev_an;
my $mean_alpha;
my $stddev_alpha;
my $mean_fpos;
my $stddev_fpos;
my $mean_hill;
my $stddev_hill;
my $mean_log_p;
my $stddev_log_p;
my $mean_r;
my $stddev_r;
my $mean_beta;
my $stddev_beta;
my $mean_betae;
my $stddev_betae;
my $mean_latent_inf;
my $stddev_latent_inf;
my $mean_exp_days;
my $stddev_exp_days;
my @pieces;

open (INFILE, "<$infile");
$line=<INFILE>;
while( $line=<INFILE>) {
    chomp($line);              # remove the newline from $line.

    @pieces = split(/ /,$line);

    if ($pieces[0] eq "an_mean") {
	$mean_an=$pieces[1];

    } elsif ($pieces[0] eq "an_std") {
	$stddev_an=$pieces[1];

    } elsif ($pieces[0] eq "alpha_mean") {
	$mean_alpha=$pieces[1];

    } elsif ($pieces[0] eq "alpha_std") {
	$stddev_alpha=$pieces[1];

    } elsif ($pieces[0] eq "beta_mean") {
	$mean_beta=$pieces[1];

    } elsif ($pieces[0] eq "beta_std") {
	$stddev_beta=$pieces[1];

    } elsif ($pieces[0] eq "betae_mean") {
	$mean_betae=$pieces[1];

    } elsif ($pieces[0] eq "betae_std") {
	$stddev_betae=$pieces[1];

    } elsif ($pieces[0] eq "fpos_mean") {
	$mean_fpos=$pieces[1];

    } elsif ($pieces[0] eq "fpos_std") {
	$stddev_fpos=$pieces[1];

    } elsif ($pieces[0] eq "kappa_mean") {
	$mean_kappa=$pieces[1];

    } elsif ($pieces[0] eq "kappa_std") {
	$stddev_kappa=$pieces[1];

    } elsif ($pieces[0] eq "hill_mean") {
	$mean_hill=$pieces[1];

    } elsif ($pieces[0] eq "hill_std") {
	$stddev_hill=$pieces[1];

    } elsif ($pieces[0] eq "r_mean") {
	$mean_r=$pieces[1];

    } elsif ($pieces[0] eq "r_std") {
	$stddev_r=$pieces[1];

    } elsif ($pieces[0] eq "log_p_mean") {
	$mean_log_p=$pieces[1];

    } elsif ($pieces[0] eq "log_p_std") {
	$stddev_log_p=$pieces[1];

    } elsif ($pieces[0] eq "exp_days_mean") {
	$mean_exp_days=$pieces[1];

    } elsif ($pieces[0] eq "exp_days_std") {
	$stddev_exp_days=$pieces[1];

    } elsif ($pieces[0] eq "latent_inf_mean") {
	$mean_latent_inf=$pieces[1];

    } elsif ($pieces[0] eq "latent_inf_std") {
	$stddev_latent_inf=$pieces[1];
    }
}
close(INFILE);



while ( $count > 0 ){
    print "$count runs remaining\n";

    open (TSTFILE, ">>$infile");
    printf TSTFILE ("PDF_on 0\n");
    my $beta= randn($mean_beta,$stddev_beta);
    if ( $model == 1 ) {
	printf TSTFILE ("infect_by_virus 0\n");
	printf TSTFILE ("beta_init %g\n",$beta);
    }
    if ( $model > 1 ) {
	printf TSTFILE ("infect_by_virus 1\n");
	printf TSTFILE ("beta_init %e\n",$beta);
    }
    if ( $model == 3 ){
	my $betae= randn($mean_betae,$stddev_betae);
	printf TSTFILE ("betae_init %e\n",$betae);
    }
    my $fpos= randn($mean_fpos,$stddev_fpos);
    printf TSTFILE ("fpos %g\n",$fpos);
    my $an= randn($mean_an,$stddev_an);
    printf TSTFILE ("an %g\n",$an);

    # r and alpha are fixed
    #my $alpha= randn($mean_alpha,$stddev_alpha);
    #printf TSTFILE ("alpha_init %g\n",$alpha);
    #my $r= randn($mean_r,$stddev_r);
    #printf TSTFILE ("r_init %g\n",$r);

    #uncomment if fitting with density dependent killing
    #printf TSTFILE ("density_killing 1\n");
    #my $hill= randn($mean_hill,$stddev_hill);
    #printf TSTFILE ("hill %g\n",$hill);

    #uncomment if fitting with T cell exhaustion term
    #my $kappa= randn($mean_kappa,$stddev_kappa);
    #printf TSTFILE ("kappa %g\n",$kappa);
    printf TSTFILE ("kappa_init 0\n");

    my $log_p= randn($mean_log_p,$stddev_log_p);
    printf TSTFILE ("log_p_init %g\n",$log_p);
    my $latent_inf= randn($mean_latent_inf,$stddev_latent_inf);
    printf TSTFILE ("latent_inf_init %g\n",$latent_inf);
    close(TSTFILE);
    system("../hhv8_sim -b -r -f $infile -c $critfile -w 1");
    $count--;
}
