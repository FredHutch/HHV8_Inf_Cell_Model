#!/usr/bin/perl
my $runs=0;
my $has_epis=1;
my @scores;
my @score1s;
my @score2s;
my @score3s;
my @score4s;
my @score5s;
my @score8s;
my @AIC_ks;
my @AIC_scores;

my @abs_scores;
my @abs_score1s;
my @abs_score2s;
my @abs_score3s;

my @shedding_scores;
my @med_hhv8_scores;
my @med_log_hhv8;
my @var_log_hhv8;
my @peak_log_hhv8;
my @init_T0;

my $crit1bin=1;
my $crit1bins=6;
my $crit2bin=0;
my $crit2bins=5;
my $crit3bin=0;
my $crit3bins=6;

my @ibirths;
my @tbirths;
my @all_exp_days;
my @alphas;
my @betas;
my @betaes;
my @log_ps;
my @rhos;
my @fposs;
my @hills;
my @ans;
my @deltas;
my @thetas;
my @kappas;
my @cd8_ic50s;
my @latent_infs;
my @rs;
my @cs;
my @infs;
my @rinfs;
my @eclipses;
my @beta_ps;
my @regs;

my $score;
my $score1;
my $score2;
my $score3;
my $alpha;
my $exp_days;
my $beta;
my $betae;
my $beta_p;
my $log_p;
my $rho;
my $fpos;
my $hill;
my $an;
my $delta;
my $theta;
my $kappa;
my $cd8_ic50;
my $latent_inf;
my $r;
my $c;
my $inf;
my $rinf;
my $eclipse;
my $reg;

my $AIC_k;

my @cat_means;
my @cat_lows;
my @cat_highs;
my @cat_acts;
my @cat_sims;

$abs_scores[$runs] = 0; 
$abs_score1s[$runs] = 0; 
$abs_score2s[$runs] = 0; 
$abs_score3s[$runs] = 0; 

$score1s[$runs] = 0; 
$score2s[$runs] = 0; 
$score3s[$runs] = 0; 
$scores[$runs] = 0; 

$shedding_scores[$runs] = 0;
$med_hhv8_scores[$runs] = 0;
$var_hhv8_scores[$runs] = 0;
$med_log_hhv8[$runs] = 0;
$var_log_hhv8[$runs] = 0;
$peak_log_hhv8[$runs] = 0;
$dur_mean_svar_cores[$runs] = 0;

my $txtfile = $ARGV[0];
my $critfile = $ARGV[1];
my $critfile2 = $ARGV[2];
my $ptid = $ARGV[3];

#70 72 71.9706
#8 9 8.89433
#10 12 11.1383
#3 6 4.77356
#1 3 2.36638
#0.5 1.5 0.856793
#40 45 42.246
#30 35 33.1551
#15 20 17.6471
#4 6 4.81283
#1.5 2 2.13904
#23.6842	28.9474	26.3158
#11.8421	14.4737	13.1579
#11.8421	14.4737	13.1579
#18.9474	23.1579	21.0526
#18.9474	23.1579	21.0526
#04.73684	05.78947	05.26316
my @pieces;
open(CRIT,"<$critfile") || die " Could not open criteria file $critfile!\n";
my $line;
my $category = 1;
my $crit_bins = $crit1bins;
while( $category <= 3) {

    my $tot_means=0;
    for (my $i=0; $i < $crit_bins; $i++) {
	$line =<CRIT>;
	chomp($line);              # remove the newline from $line.

	@pieces = split(/\t/,$line);
	$cat_lows[$category-1][$i]=$pieces[0];
	$cat_highs[$category-1][$i]=$pieces[1];
	$cat_acts[$category-1][$i]=$pieces[2];
	my $mean = $pieces[2];

	# skip 1st entry (percent neg swabs)
	if ($category != 1 || $i > 0) {
	    $tot_means+=$mean;
	    my $index = $i + 1;
	    print STDERR "Category $category [$i]: low=$cat_lows[$category-1][$i], high=$cat_highs[$category-1][$i], act=$cat_acts[$category-1][$i]\n";
	}
    }
    if ($crit_bins > 0) {
	if ($category == 1) {
	    $cat_means[$category-1]=$tot_means/($crit_bins-1);
	} else {
	    $cat_means[$category-1]=$tot_means/$crit_bins;
	}
	print STDERR "Category $category: $crit_bins bins, mean=$cat_means[$category-1]\n";
    } else {
	$cat_means[$category-1]=0;
    }

    $category++;
    if ($category == 2) {
	$crit_bins = $crit2bins;
    } else {
	$crit_bins = $crit3bins;
    }
}
close(CRIT);
open(CRIT2,"<$critfile2") || die " Could not open criteria file $critfile2!\n";

# read swab percentage, median and peak log hhv8
$line=<CRIT2>;
chomp($line);              # remove the newline from $line.
@pieces = split(/[\t ]/,$line);
my $shedding=$pieces[0];
my $med_hhv8_target=$pieces[1];
my $var_hhv8_target=$pieces[2];
my $peak_hhv8_target=$pieces[3];
printf STDERR ("Shedding target = %g\n",$shedding);
printf STDERR ("Median HHV8 target = %g\n",$med_hhv8_target);
printf STDERR ("Variance HHV8 target = %g\n",$var_hhv8_target);
printf STDERR ("Peak HHV8 target = %g\n",$peak_hhv8_target);
close(CRIT2);

open(TXT,"<$txtfile") || die " Could not open input file $txtfile!\n";

$line=<TXT>;
while( $line=<TXT>) {
    chomp($line);              # remove the newline from $line.

    @pieces = split(/ /,$line);


    if ($pieces[1] eq "Regions") {
	$reg=$pieces[3];

    } elsif ($pieces[1] eq "N_runs") {
	$N_runs=$pieces[3];

    } elsif ($pieces[1] eq "beta_init") {
	$beta=$pieces[3];

    } elsif ($pieces[1] eq "kappa_init") {
	$kappa=$pieces[3];

    } elsif ($pieces[1] eq "cd8_ic50_init") {
	$cd8_ic50=$pieces[3];

    } elsif ($pieces[1] eq "log_p_init") {
	$log_p=$pieces[3];

    } elsif ($pieces[1] eq "c_init") {
	$c=$pieces[3];

    } elsif ($pieces[1] eq "delta_init") {
	$delta=$pieces[3];

    } elsif ($pieces[1] eq "r_init") {
	$r=$pieces[3];

    } elsif ($pieces[1] eq "inf_init") {
	$inf=$pieces[3];

    } elsif ($pieces[1] eq "rinf_init") {
	$rinf=$pieces[3];

    } elsif ($pieces[1] eq "eclipse_init") {
	$eclipse=$pieces[3];

    } elsif ($pieces[1] eq "rho_init") {
	$rho=$pieces[3];

    } elsif ($pieces[1] eq "an") {
	$an=$pieces[3];

    } elsif ($pieces[1] eq "alpha_init") {
	$alpha=$pieces[3];

    } elsif ($pieces[1] eq "betae_init") {
	$betae=$pieces[3];

    } elsif ($pieces[1] eq "exp_days_init") {
	$exp_days=$pieces[3];

    } elsif ($pieces[1] eq "AIC_K") {
	$AIC_k = $pieces[3];

# Initial T0 is 2206 (inf reg=47968 (3 at end), hhv8 reg=85620 (7 at end)). Beginning simulation!
    } elsif ($pieces[0] eq "Initial" && $pieces[1] eq "T0") {
	$init_T0[$runs] += $pieces[3];

#Total infect cell births 26
    } elsif ($pieces[0] eq "Total" && $pieces[3] eq "births") {
	$ibirths[$runs] += $pieces[4];
 
#Total T cells produced 1259
    } elsif ($pieces[0] eq "Total" && $pieces[1] eq "T" && $pieces[3] eq "produced") {
	$tbirths[$runs] += $pieces[4];
 
#Episode end at 178.000000 days: Ve peak=1047263 (bin=4), Episode duration=28.000000 days (peak at t=152.050000)
    } elsif ($pieces[0] eq "Episode" && $pieces[1] eq "end") {
	@subpieces = split(/=/,$pieces[6]);
	my $peak = log($subpieces[1])/log(10);
	#if ($peak > $peak_log_hhv8[$runs]) {
	#    $peak_log_hhv8[$runs] = $peak;
	#}
    } elsif ($pieces[0] eq "Median" && $pieces[2] eq "swab") {
	$med_log_hhv8[$runs] = $pieces[5];
	if ($pieces[5] > 0 || $med_hhv8_target > 0) {
	    $med_hhv8_scores[$runs] = abs($pieces[5] - $med_hhv8_target)/($med_hhv8_target);
	} else {
	    $med_hhv8_scores[$runs] = 0;
	}
	printf STDERR ("For set %d: med_hhv8_score = %g (%g vs %g)\n",
		$runs+1,$med_hhv8_scores[$runs],$pieces[5],$med_hhv8_target);
    
    } elsif ($pieces[0] eq "Mean" && $pieces[2] eq "swab") {
	printf STDERR ("For set %d: swab_mean = %g\n",
		$runs+1,$pieces[5]);
    
    } elsif ($pieces[0] eq "Var" && $pieces[2] eq "swab") {
	$var_log_hhv8[$runs] = $pieces[5];
	if ($pieces[5] > 0 || $var_hhv8_target > 0) {
	    $var_hhv8_scores[$runs] = 2*abs($pieces[5] - $var_hhv8_target)/($var_hhv8_target+$pieces[5]);
	} else {
	    $var_hhv8_scores[$runs] = 0;
	}
	printf STDERR ("For set %d: var_hhv8_score = %g (%g vs %g)\n",
		$runs+1,$var_hhv8_scores[$runs],$pieces[5],$var_hhv8_target);
    
    } elsif ($pieces[0] eq "Mean" && $pieces[2] eq "peak") {
	printf STDERR ("For set %d: peak_mean = %g\n",
		$runs+1,$pieces[5]);
	$peak_log_hhv8[$runs] = $pieces[5];
	if ($peak_log_hhv8[$runs] > 0) {
	    $peak_hhv8_scores[$runs] = abs($pieces[5] - $peak_hhv8_target)/($peak_hhv8_target);
	} else {
	    $peak_hhv8_scores[$runs] = 0;
	}
	printf STDERR ("For set %d: peak_hhv8_score = %g (%g vs %g)\n",
		$runs+1,$peak_hhv8_scores[$runs],$peak_log_hhv8[$runs],$peak_hhv8_target);
    
    } elsif ($pieces[0] eq "Var" && $pieces[1] eq "peak") {
	printf STDERR ("For set %d: peak_var = %g\n",
		$runs+1,$pieces[4]);

    } elsif ($pieces[0] eq "Mean" && $pieces[2] eq "duration") {
	printf STDERR ("For set %d: dur_mean = %g\n",
		$runs+1,$pieces[4]);
    
    } elsif ($pieces[0] eq "Var" && $pieces[2] eq "duration") {
	printf STDERR ("For set %d: dur_var = %g\n",
		$runs+1,$pieces[4]);

    } elsif ($pieces[0] eq "No" && $pieces[1] eq "episodes!") {
	$has_epis=0;

#Mean log swab VL = 5.299265
#Var log swab VL = 0.230047s
#Mean log peak VL = 5.447016
#Var peak VL = 0.173452
#Mean episode duration = 0.368421 days
#Var episode duration = 10.180055 days
#Final percentage of swabs that were positive= 13.282407
    } elsif ($pieces[0] eq "Final") {
	$shedrates[$runs]=$pieces[7];
	if ($pieces[7] > 0 || $shedding > 0) {
	    $shedding_scores[$runs] = abs($pieces[7] - $shedding)/($shedding);
	} else {
	    $shedding_scores[$runs] = 0;
	}
	printf STDERR ("For set %d: shed_score = %g (%g vs %g)\n",
		$runs+1,$shedding_scores[$runs],$pieces[7],$shedding);
	$init_T0[$runs] = $init_T0[$runs] / $N_runs;
	$ibirths[$runs] = $ibirths[$runs] / $N_runs;
	$tbirths[$runs] = $tbirths[$runs] / $N_runs;
	if ($has_epis == 1) {
	    $runs++;
	}
	$init_T0[$runs] = 0;
	$ibirths[$runs] = 0;
	$tbirths[$runs] = 0;
	$peak_log_hhv8[$runs] = 0;

	$has_epis = 1;
	$score1s[$runs] = 0; 
	$score2s[$runs] = 0; 
	$score3s[$runs] = 0; 
	$scores[$runs] = 0; 

	$abs_scores[$runs] = 0; 
	$abs_score1s[$runs] = 0; 
	$abs_score2s[$runs] = 0; 
	$abs_score3s[$runs] = 0; 

	$crit1bin=1;
	$crit2bin=0;
	$crit3bin=0;
	$shedding_scores[$runs] = 0;
    } elsif ($pieces[0] eq "Score1") {
	$score1 = $pieces[2];
	$score1s[$runs] = $score1; 

    } elsif ($pieces[0] eq "Score2") {
	$score2 = $pieces[2];
	$score2s[$runs] = $score2; 

    } elsif ($pieces[0] eq "Score3") {
	$score3 = $pieces[2];
	$score3s[$runs] = $score3; 

#criteria 1[2]: 4.824561 vs. mean of 8.894330 - err = 4.069769 * 0.029731 * 1.000000
    } elsif ($pieces[0] eq "criteria" && $pieces[1] =~ /^1/) {
#
	    $crit1[$runs][$crit1bin]=$pieces[2];
	    $sim1_vals[$runs][$crit1bin]=$pieces[2];
	    $abs_score1s[$runs] += abs($pieces[2] - $cat_acts[0][$crit1bin])/($cat_means[0]*($crit1bins-1));
	    #print STDERR "Category 1: [$crit1bin]: $pieces[2] vs $cat_acts[0][$crit1bin] err=$abs_score1s[$runs]\n";
	    $crit1bin++;

#criteria 2[1]: 10.714286 vs. mean of 42.246000 - err = 31.531714 * 1.000000
    } elsif ($pieces[0] eq "criteria" && $pieces[1] =~ /^2/) {
	    $crit2[$runs][$crit2bin]=$pieces[2];
	    $sim2_vals[$runs][$crit2bin]=$pieces[2];
	    $abs_score2s[$runs] += abs($pieces[2] - $cat_acts[1][$crit2bin])/($cat_means[1]*$crit2bins);
	    #print STDERR "Category 2: [$crit2bin]: $pieces[2] vs $cat_acts[1][$crit2bin] err=$abs_score2s[$runs]\n";
	    $crit2bin++;

#criteria 3[1]: 35.526316 vs. mean of 26.315800 - err = 9.210516 * 1.000000
    } elsif ($pieces[0] eq "criteria" && $pieces[1] =~ /^3/) {
	    $crit3[$runs][$crit3bin]=$pieces[2];
	    $sim3_vals[$runs][$crit3bin]=$pieces[2];
	    $abs_score3s[$runs] += abs($pieces[2] - $cat_acts[2][$crit3bin])/($cat_means[2]*$crit3bins);
	    #print STDERR "Category 3: [$crit3bin]: $pieces[2] vs $cat_acts[2][$crit3bin] err=$abs_score3s[$runs]\n";
	    $crit3bin++;

#Total Score = 0.450561 beta=153.317000e-8 latent_inf=0.082958 log_p=4.035470 c=10.022200 theta=4.554530 delta=0.000081 r=36.095500 inf=122358.000000 rinf=362272.000000 beta_e=88.511510e-11 rho=0.664769 eclipse=0.828496 fpos=0.490412 an=0.000000 time=2128.000000
#Total Score = 0.534799 beta=297859475.909418e-8 latent_inf=0.832041 log_p=4.000000 c=6.000000 theta=5.000000 delta=0.001000 r=500.000000 inf=0.000000 rinf=116685.000000 rho=0.763742 eclipse=0.850000 fpos=1.316655 an=1.500000 cd8_ic50=10000.000000 kappa=5.000000 alpha=0.01 exp_days=0 hill=0.01 time=700.000000
    } elsif ($pieces[1] eq "Score") {
	$score = $pieces[3];

	$scores[$runs]=$score;

	#@subpieces = split(/=/,$pieces[4]);
	#$beta = $subpieces[1];

	@subpieces = split(/=/,$pieces[5]);
	$latent_inf = $subpieces[1];

	@subpieces = split(/=/,$pieces[8]);
	$theta = $subpieces[1];

	@subpieces = split(/=/,$pieces[15]);
	$fpos = $subpieces[1];

	@subpieces = split(/=/,$pieces[21]);
	$hill = $subpieces[1];

	$AIC_ks[$runs] = $AIC_k;
	$AIC_scores[$runs] = ($abs_score1s[$runs]+$abs_score2s[$runs]+$abs_score3s[$runs]);
	$AIC_scores[$runs] += 2*$AIC_ks[$runs];

	$abs_scores[$runs]=1.0 / ($abs_score1s[$runs]+$abs_score2s[$runs]+$abs_score3s[$runs]);

	$beta_p = $beta * (10 ** $log_p);

	$betas[$runs]=$beta;
	$rhos[$runs]=$rho;
	$fposs[$runs]=$fpos;
	$hills[$runs]=$hill;
	$ans[$runs]=$an;
	$infs[$runs]=$inf;
	$rinfs[$runs]=$rinf;
	$eclipses[$runs]=$eclipse;
	$log_ps[$runs]=$log_p;
	$rs[$runs]=$r;
	$cs[$runs]=$c;
	$thetas[$runs]=$theta;
	$kappas[$runs]=$kappa;
	$cd8_ic50s[$runs]=$cd8_ic50;
	$deltas[$runs]=$delta;
	$latent_infs[$runs]=$latent_inf;
	$beta_ps[$runs]=$beta_p;
	$regs[$runs]=$reg;
	$alphas[$runs]=$alpha;
	$betaes[$runs]=$betae;
	$all_exp_days[$runs]=$exp_days;
    }
}
close(TXT);

printf ("Run,mean beta,mean log_p,mean latent_inf,mean theta,mean delta,mean r,mean c,mean rho,mean inf,mean rinf,mean eclipse,mean hill,mean To,mean fpos,mean an,regs");

printf (",avg shedrate,score1,score2,score3,score");
for (my $j=1; $j < $crit1bins;$j++){
    printf (",sim1_%d,act1_%d",$j,$j);
}
printf (",new_score1");
for (my $j=0; $j < $crit2bins;$j++){
    printf (",sim2_%d,act2_%d",$j,$j);
}
printf (",new_score2");
for (my $j=0; $j < $crit3bins;$j++){
    printf (",sim3_%d,act3_%d",$j,$j);
}
printf (",new_score3,new_score,AIC_k,AIC_score,mean kappa,mean cd8_ic50,total new score");
printf (",old score,shedding,shedding target,shedding score");
printf (",med log hhv8,med log hhv8 target,med hhv8 score");
printf (",var log hhv8, var target,var hhv8 score,ptid,exp days,alpha,avg ibirths,avg tbirths");
printf (",peak log hhv8, peak target,peak hhv8 score,three score,betae\n");
for (my $i=0; $i < $runs; $i++){
    printf ("%d,%e,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%e,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf",
	$i+1,$betas[$i],$log_ps[$i],$latent_infs[$i],$thetas[$i],$deltas[$i],$rs[$i],$cs[$i],$rhos[$i],
	$infs[$i],$rinfs[$i],$eclipses[$i],$hills[$i],$init_T0[$i],$fposs[$i],$ans[$i],$regs[$i],
	$shedrates[$i],$score1s[$i],$score2s[$i],$score3s[$i],$scores[$i]);
    for (my $j=1; $j < $crit1bins;$j++){
	printf (",%lf,%lf", $sim1_vals[$i][$j],$cat_acts[0][$j]);
    }
    printf (",%lf",$abs_score1s[$i]);
    for (my $j=0; $j < $crit2bins;$j++){
	printf (",%lf,%lf", $sim2_vals[$i][$j],$cat_acts[1][$j]);
    }
    printf (",%lf",$abs_score2s[$i]);
    for (my $j=0; $j < $crit3bins;$j++){
	printf (",%lf,%lf", $sim3_vals[$i][$j],$cat_acts[2][$j]);
    }
    printf (",%lf",$abs_score3s[$i]);
    printf (",%lf",$abs_scores[$i]);
    printf (",%d",$AIC_k[$i]);
    printf (",%lf,%lf,%lf",$AIC_scores[$i],$kappas[$i],$cd8_ic50s[$i]);
    printf (",%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%e\n",
	($shedding_scores[$i]+$med_hhv8_scores[$i])/2,
	($shedding_scores[$i]+$med_hhv8_scores[$i]+$var_hhv8_scores[$i])/3,
        $shedrates[$i],$shedding,$shedding_scores[$i],
        $med_log_hhv8[$i],$med_hhv8_target,$med_hhv8_scores[$i],
	$var_log_hhv8[$i],$var_hhv8_target,$var_hhv8_scores[$i],$ptid,
	$all_exp_days[$i],$alphas[$i],$ibirths[$i],$tbirths[$i],
        $peak_log_hhv8[$i],$peak_hhv8_target,$peak_hhv8_scores[$i],
	($shedding_scores[$i]+$med_hhv8_scores[$i]+$peak_hhv8_scores[$i])/3,$betaes[$i]);
}
