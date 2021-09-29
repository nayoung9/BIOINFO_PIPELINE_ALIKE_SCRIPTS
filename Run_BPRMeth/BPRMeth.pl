#!/usr/bin/perl
#
use strict;
use warnings;
use FindBin qw($Bin);


my $f_all_methylcalls = shift;
my $outdir = shift;
my $annotation_file = "/mss2/projects/Minipig_Methylome/Analysis_scripts/Clustering/pig_anno.bed";

#input file process
`mkdir -p $outdir`;
`cut -f1,2,4 $f_all_methylcalls > $outdir/tmpinput`;
`cut -f3 $f_all_methylcalls > $outdir/tmpinput2`;
`paste $outdir/tmpinput $outdir/tmpinput2 > $outdir/met_input.txt`;
`rm -rf $outdir/tmpinput; rm -rf $outdir/tmpinput2`;

my $methylinput_file = "$outdir/met_input.txt";
#run clustering 

`Rscript  $Bin/BPRMeth.R $annotation_file $methylinput_file $outdir &> $outdir/log`;
