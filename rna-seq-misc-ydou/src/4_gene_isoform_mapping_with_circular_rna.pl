#!/usr/bin/perl -w

use Getopt::Long;

my ($in_file, $out_file, $help);

Getopt::Long::GetOptions (
  'in|I=s'          =>  \$in_file,
  'out|O=s'         =>  \$out_file,
  'help|H!'         =>  \$help
);


if (defined($help)) {
  print "
Usage: ./4_gene_isoform_mapping_with_circular_rna.pl  -I in.gtf  -O output.txt

Arguments:

    -I, --in
          input gtf file
    -O, --out
          output transcript and gene mapping table file
\n";
	exit 0;
}

if ( !defined($in_file) or !defined($out_file) )  {
  print "Please use the --help or -H option to get usage information.\n";
  exit 1;
} 

open(IN, $in_file);
open(OUT,">", $out_file);
%mapping=();
while($raw=<IN>)
{
	$raw=~s/\n|\r|\;|\"//g;
	@inf=split(/\s+/,$raw);
	$mapping{$inf[11]}=$inf[9];

}
close(IN);
foreach $key (sort keys(%mapping))
{
	print OUT "$mapping{$key}\t$key\n";
}
close(OUT);
