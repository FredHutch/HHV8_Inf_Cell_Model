#!/bin/sh
#./indiv_rollup.pl init_arms_5.txt both_pos.crit both_pos_summary.crit 0 > both_pos_5.csv
#head -n 1 both_pos_5.csv > top100_both_pos_5.csv
#grep -v Run both_pos_5.csv > temp_scores.csv
#sort -n -k82  -t',' temp_scores.csv |head -n 100 >>  top100_both_pos_5.csv
#./top100_dist.pl top100_both_pos_5.csv 5 >> both_pos_5.in
#tail -n 20 both_pos_5.in
#cut -f82 -d',' top100_both_pos_5.csv | head -n 2
#cut -f82 -d',' top100_both_pos_5.csv | tail -n 1
echo "For both neg..."
./indiv_rollup.pl init_arms_5.txt both_neg.crit both_neg_summary.crit 0 > both_neg_5.csv
head -n 1 both_neg_5.csv > top100_both_neg_5.csv
grep -v Run both_neg_5.csv > temp_scores.csv
sort -n -k82  -t',' temp_scores.csv |head -n 100 >>  top100_both_neg_5.csv
./top100_dist.pl top100_both_neg_5.csv 5 >> both_neg_5.in
tail -n 20 both_neg_5.in
cut -f82 -d',' top100_both_neg_5.csv | head -n 2
cut -f82 -d',' top100_both_neg_5.csv | tail -n 1
echo
echo "For HIV pos..."
./indiv_rollup.pl init_arms_5.txt hiv_pos.crit hiv_pos_summary.crit 0 > hiv_pos_5.csv
head -n 1 hiv_pos_5.csv > top100_hiv_pos_5.csv
grep -v Run hiv_pos_5.csv > temp_scores.csv
sort -n -k82  -t',' temp_scores.csv |head -n 100 >>  top100_hiv_pos_5.csv
./top100_dist.pl top100_hiv_pos_5.csv 5 >> hiv_pos_5.in
tail -n 20 hiv_pos_5.in
cut -f82 -d',' top100_hiv_pos_5.csv | head -n 2
cut -f82 -d',' top100_hiv_pos_5.csv | tail -n 1
echo
echo "For KS pos..."
./indiv_rollup.pl init_arms_5.txt ks_pos.crit ks_pos_summary.crit 0 > ks_pos_5.csv
head -n 1 ks_pos_5.csv > top100_ks_pos_5.csv
grep -v Run ks_pos_5.csv > temp_scores.csv
sort -n -k82  -t',' temp_scores.csv |head -n 100 >>  top100_ks_pos_5.csv
./top100_dist.pl top100_ks_pos_5.csv 5 >> ks_pos_5.in
tail -n 20 ks_pos_5.in
cut -f82 -d',' top100_ks_pos_5.csv | head -n 2
cut -f82 -d',' top100_ks_pos_5.csv | tail -n 1
