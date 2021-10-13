#!/usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename; 
use Cwd 'abs_path';

my $motif_dir = shift;
my $f_inputfasta = shift;
my $cutoff = shift;
my $out_dir = shift;

`mkdir -p $out_dir`;

my @ar_motifs = <$motif_dir/*meme>;
foreach my $thismotif (@ar_motifs){
	if ($thismotif !~ /meme$/){next;}
	print $thismotif."\n";
	my $prefix = basename($thismotif,".meme");
	my $cmd = "fimo --qv-thresh --thresh $cutoff  --oc $out_dir/$prefix $thismotif $f_inputfasta";
	#print STDERR $cmd."\n";
	`$cmd`;
}

my $cmd_summary = "cat $out_dir/*/*.tsv | egrep -v \"#|motif\" > $out_dir.result.txt";
`$cmd_summary`;

