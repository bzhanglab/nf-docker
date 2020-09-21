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
Usage: ./3_circular_linear_to_psedo_linear.pl -I in.fa  -O output.fa

Arguments:

    -I, --in
          input transcript file
    -O, --out
          output transcript file with linear transcript of circRNA 
		  transferred to psedo linear transcript
\n";
	exit 0;
}

if ( !defined($in_file) or !defined($out_file) )  {
  print "Please use the --help or -H option to get usage information.\n";
  exit 1;
} 


$read_length=75;
open(IN, $in_file);
open(OUT,">", $out_file);

while($raw=<IN>)
{
	$raw2=<IN>;
	$raw=~s/\n|\r//g;
	$raw2=~s/\n|\r//g;
	if($raw=~/circ_/)
	{
		if(length($raw2)>$read_length&&$raw2!~/N/)
		{
			print OUT $raw,"\n";
			$head_str=substr($raw2,0,75);
			print OUT "$raw2$head_str\n";
		}
	}
	else
	{
		if(length($raw2)>$read_length&&$raw2!~/N/)
		{
			print OUT $raw,"\n";
			print OUT $raw2,"\n";
		}
	}
}
close(IN);
close(OUT);
