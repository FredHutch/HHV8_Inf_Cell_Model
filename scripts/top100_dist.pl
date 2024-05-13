#!/usr/bin/perl
my $txtfile = $ARGV[0];
my $model = $ARGV[1];

# use sorted file and select top 100 entries to get distributions
open(TXT,"<$txtfile") || die " Could not open input file $txtfile!\n";

#Run,mean beta,mean log_p,mean latent_inf,mean theta,mean delta,mean r,mean c,mean rho,mean inf,mean rinf,mean eclipse,mean hill,mean To,mean fpos,mean an,regs,avg shedrate,score1,score2,score3,score,sim1_1,act1_1,sim1_2,act1_2,sim1_3,act1_3,sim1_4,act1_4,sim1_5,act1_5,new_score1,sim2_0,act2_0,sim2_1,act2_1,sim2_2,act2_2,sim2_3,act2_3,sim2_4,act2_4,new_score2,sim3_0,act3_0,sim3_1,act3_1,sim3_2,act3_2,sim3_3,act3_3,sim3_4,act3_4,sim3_5,act3_5,new_score3,new_score,AIC_k,AIC_score,mean kappa,mean cd8_ic50,total new score,old score,shedding,shedding target,shedding score,med log hhv8,med log hhv8 target,med hhv8 score,var log hhv8, var target,var hhv8 score,ptid,exp days,alpha,avg ibirths,avg tbirths,peak log hhv8, peak target,peak hhv8 score,three score
my $num_scores=0;
my $sum_scores=0;
my @all_scores;
my @raw_scores;
my @alpha;
my @beta;
my @betae;
my @latent_inf;
my @hill;
my @fpos;
my @log_p;
my @exp_days;
my @r;
my @an;
my @shedding;
my @ptid;

$line=<TXT>;	# skip header

while( $line=<TXT>) {
    chomp($line);              # remove the newline from $line.
    my @pieces = split(/,/,$line);

    if ($num_scores < 100) {
	$raw_scores[$num_scores]=$pieces[81];
	$all_scores[$num_scores]=$pieces[81];
	$latent_inf[$num_scores]=$pieces[3];
	$beta[$num_scores]=$pieces[1];
	$betae[$num_scores]=$pieces[82];
	$log_p[$num_scores]=$pieces[2];
	$an[$num_scores]=$pieces[15];
	$hill[$num_scores]=$pieces[12];
	$fpos[$num_scores]=$pieces[14];
	$r[$num_scores]=$pieces[6];
	$ptid[$num_scores]=$pieces[73];
	$exp_days[$num_scores]=$pieces[74];
	$alpha[$num_scores]=$pieces[75];
	$sum_scores+=$all_scores[$num_scores];
	$shedding[$num_scores]=$pieces[17];
	#printf STDERR ("Score %d = %g\n",$num_scores+1,$all_scores[$num_scores]);
	$num_scores++;
    }
}
close(TXT);

my $summed_beta=0;
my $summed_betae=0;
my $summed_an=0;
my $summed_hill=0;
my $summed_fpos=0;
my $summed_r=0;
my $summed_log_p=0;
my $summed_latent_inf=0;
my $summed_exp_days=0;
my $summed_alpha=0;
my $low_shedders=0;

for (my $i=0; $i < $num_scores; $i++) {
    #printf STDERR ("Pick %d (val=%g, beta=%g,fpos=%g,lat_inf=%g)\n",$i,$all_scores[$i],$beta[$i],$fpos[$i],$latent_inf[$i]);
    $summed_beta+=$beta[$i];
    $summed_betae+=$betae[$i];
    $summed_an+=$an[$i];
    $summed_hill+=$hill[$i];
    $summed_fpos+=$fpos[$i];
    $summed_r+=$r[$i];
    $summed_log_p+=$log_p[$i];
    $summed_latent_inf+=$latent_inf[$i];
    $summed_exp_days+=$exp_days[$i];
    $summed_alpha+=$alpha[$i];
    if ($shedding[$i] < 55) {
	$low_shedders++;
	$low_scores += $raw_scores[$i];
    } else {
	$high_shedders++;
	$high_scores += $raw_scores[$i];
    }
}
my $beta_m= $summed_beta/$num_scores;
my $beta_sqtotal = 0;
my $betae_m= $summed_betae/$num_scores;
my $betae_sqtotal = 0;
my $an_m= $summed_an/$num_scores;
my $an_sqtotal = 0;
my $hill_m= $summed_hill/$num_scores;
my $fpos_m= $summed_fpos/$num_scores;
my $r_m= $summed_r/$num_scores;
my $log_p_m= $summed_log_p/$num_scores;
my $hill_sqtotal = 0;
my $fpos_sqtotal = 0;
my $log_p_sqtotal = 0;
my $r_sqtotal = 0;
my $latent_inf_m= $summed_latent_inf/$num_scores;
my $latent_inf_sqtotal = 0;
my $exp_days_m= $summed_exp_days/$num_scores;
my $exp_days_sqtotal = 0;
my $alpha_m= $summed_alpha/$num_scores;
my $alpha_sqtotal = 0;
for (my $i=0; $i < $num_scores; $i++) {
    $beta_sqtotal += ($beta_m-$beta[$i]) ** 2;
    $betae_sqtotal += ($betae_m-$betae[$i]) ** 2;
    $an_sqtotal += ($an_m-$an[$i]) ** 2;
    $hill_sqtotal += ($hill_m-$hill[$i]) ** 2;
    $fpos_sqtotal += ($fpos_m-$fpos[$i]) ** 2;
    $r_sqtotal += ($r_m-$r[$i]) ** 2;
    $log_p_sqtotal += ($log_p_m-$log_p[$i]) ** 2;
    $latent_inf_sqtotal += ($latent_inf_m-$latent_inf[$i]) ** 2;
    $exp_days_sqtotal += ($exp_days_m-$exp_days[$i]) ** 2;
    $alpha_sqtotal += ($alpha_m-$alpha[$i]) ** 2;
}
my $beta_s= ($beta_sqtotal / ($num_scores-1)) ** 0.5;
my $betae_s= ($betae_sqtotal / ($num_scores-1)) ** 0.5;
my $an_s= ($an_sqtotal / ($num_scores-1)) ** 0.5;
my $log_p_s= ($log_p_sqtotal / ($num_scores-1)) ** 0.5;
my $an_s= ($an_sqtotal / ($num_scores-1)) ** 0.5;
my $log_p_s= ($log_p_sqtotal / ($num_scores-1)) ** 0.5;
my $hill_s= ($hill_sqtotal / ($num_scores-1)) ** 0.5;
my $fpos_s= ($fpos_sqtotal / ($num_scores-1)) ** 0.5;
my $r_s= ($r_sqtotal / ($num_scores-1)) ** 0.5;
my $latent_inf_s= ($latent_inf_sqtotal / ($num_scores-1)) ** 0.5;
my $exp_days_s= ($exp_days_sqtotal / ($num_scores-1)) ** 0.5;
my $alpha_s= ($alpha_sqtotal / ($num_scores-1)) ** 0.5;
my @sorted_betas = sort { $a <=> $b } @beta;
my @sorted_betaes = sort { $a <=> $b } @betae;
my @sorted_an = sort { $a <=> $b } @an;
my @sorted_fpos = sort { $a <=> $b } @fpos;
my @sorted_an = sort { $a <=> $b } @an;
my @sorted_hill = sort { $a <=> $b } @hill;
my @sorted_fpos = sort { $a <=> $b } @fpos;
my @sorted_log_p = sort { $a <=> $b } @log_p;
my @sorted_r = sort { $a <=> $b } @r;
my @sorted_latent_infs = sort { $a <=> $b } @latent_inf;
my @sorted_exp_days = sort { $a <=> $b } @exp_days;
my @sorted_alpha = sort { $a <=> $b } @alpha;
printf STDERR ("expand_days: median = %g, mean = %g, stddev = %g\n",$sorted_exp_days[($num_scores/2)],$exp_days_m,$exp_days_s);
printf STDERR ("alpha: median = %g, mean = %g, stddev = %g\n",$sorted_alpha[($num_scores/2)],$alpha_m,$alpha_s);
printf STDERR ("beta: median = %g, mean = %g, stddev = %g\n",$sorted_betas[($num_scores/2)],$beta_m,$beta_s);
printf STDERR ("betae: median = %g, mean = %g, stddev = %g\n",$sorted_betaes[($num_scores/2)],$betae_m,$betae_s);
printf STDERR ("an: median = %g, mean = %g, stddev = %g\n",$sorted_an[($num_scores/2)],$an_m,$an_s);
printf STDERR ("hill: median = %g, mean = %g, stddev = %g\n",$sorted_hill[($num_scores/2)],$hill_m,$hill_s);
printf STDERR ("fpos: median = %g, mean = %g, stddev = %g\n",$sorted_fpos[($num_scores/2)],$fpos_m,$fpos_s);
printf STDERR ("r: median = %g, mean = %g, stddev = %g\n",$sorted_r[($num_scores/2)],$r_m,$r_s);
printf STDERR ("log_p: median = %g, mean = %g, stddev = %g\n",$sorted_log_p[($num_scores/2)],$log_p_m,$log_p_s);
printf STDERR ("latent_inf: median = %g, mean = %g, stddev = %g\n",$sorted_latent_infs[($num_scores/2)],$latent_inf_m,$latent_inf_s);
if ($low_shedders > 0) {
    printf STDERR ("Percent low shedders = %g (avg score=%g)\n",100*$low_shedders/$num_scores,$low_scores/$low_shedders);
} else {
    printf STDERR ("No low shedders\n");
}
if ($high_shedders > 0) {
    printf STDERR ("Percent high shedders = %g (avg score=%g)\n",100*$high_shedders/$num_scores,$high_scores/$high_shedders);
} else {
    printf STDERR ("No high shedders\n");
}
printf ("PDF_on 2\n");
if ( $model > 8) {
    printf ("beta_mean %e\n",$beta_m);
    printf ("beta_std %e\n",$beta_s);
}
if ( $model == 2 || $model == 4 || $model == 5 || $model == 8) {
    printf ("beta_mean %g\n",$beta_m);
    printf ("beta_std %g\n",$beta_s);
}
if ( $model == 10 ) {
    printf ("betae_mean %e\n",$betae_m);
    printf ("betae_std %e\n",$betae_s);
}
if ( $model == 1 || $model == 2 || $model >= 5) {
    printf ("fpos_mean %g\n",$fpos_m);
    printf ("fpos_std %g\n",$fpos_s);
}
if ( $model == 3 || $model == 4 || $model == 5 || $model >= 8) {
    printf ("an_mean %g\n",$an_m);
    printf ("an_std %g\n",$an_s);
}
if ( $model == 7) {
    printf ("hill_mean %g\n",$hill_m);
    printf ("hill_std %g\n",$hill_s);
}
#printf ("exp_days_mean %g\n",$exp_days_m);
#printf ("exp_days_std %g\n",$exp_days_s);
if ( $model <= 5) {
    printf ("alpha_mean %g\n",$alpha_m);
    printf ("alpha_std %g\n",$alpha_s);
    printf ("r_mean %g\n",$r_m);
    printf ("r_std %g\n",$r_s);
}
if ( $model != 6) {
    printf ("kappa_mean 0\n");
}
printf ("log_p_mean %g\n",$log_p_m);
printf ("log_p_std %g\n",$log_p_s);
printf ("latent_inf_mean %g\n",$latent_inf_m);
printf ("latent_inf_std %g\n",$latent_inf_s);
