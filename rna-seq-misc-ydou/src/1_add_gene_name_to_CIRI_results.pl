#!/usr/bin/perl -w
use Getopt::Long;

my ($ref_file, $in_file, $out_file, $help);

Getopt::Long::GetOptions (
  'in|I=s'          =>  \$in_file,
  'ref_file|F=s'        =>  \$ref_file,
  'out|O=s'         =>  \$out_file,
  'help|H!'         =>  \$help
);


if (defined($help)) {
  print "
Usage: ./1_add_gene_name_to_CIRI_results.pl -I in.txt -F ref.fa -O output.txt

Arguments:

    -I, --in
          input CIRI results
    -O, --out
          output CIRI results with gene names added
    -F, --ref_file 
	      reference gtf file

\n";
	exit 0;
}

if ( !defined($in_file) or !defined($out_file) or !defined($ref_file) ) {
  print "Please use the --help or -H option to get usage information.\n";
  exit 1;
} 

%id_to_name=();


open(IN, $ref_file);
while($raw=<IN>)
{
	$raw=~s/\n|\r|"|\;//g;
	@inf=split(/\s+/,$raw);
	$id_to_name{$inf[11]}=$inf[9];
}
close(IN);

open(IN, $in_file);
open(OUT, ">", $out_file);
$raw=<IN>;
print OUT $raw;
while($raw=<IN>)
{
	$raw=~s/\n|\r//g;
	@inf=split(/\t/,$raw);
	if($inf[9]!~/_/)
	{
		print OUT $raw,"\n";
	}
	else
	{
		@inf2=split(/,/,$inf[9]);
		print OUT "$inf[0]\t$inf[1]\t$inf[2]\t$inf[3]\t$inf[4]\t$inf[5]\t$inf[6]\t$inf[7]\t$inf[8]";
		print OUT "\t$id_to_name{$inf2[0]}";
		print OUT "\t$inf[10]\t$inf[11]\n";
	}
}
close(IN);
close(OUT);
