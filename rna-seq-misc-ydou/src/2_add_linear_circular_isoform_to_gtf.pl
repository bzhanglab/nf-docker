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
Usage: ./2_add_linear_circular_isoform_to_gtf.pl -I CIRI_results.txt -F in_ref.gtf -O out_ref.gtf

Arguments:

    -I, --in
          input CIRI results with added gene names
    -O, --out
          output reference gtf with circular RNA
    -F, --ref_file 
	      reference gtf file

\n";
	exit 0;
}

if ( !defined($in_file) or !defined($out_file) or !defined($ref_file) ) {
  print "Please use the --help or -H option to get usage information.\n";
  exit 1;
} 


$supporting_read=10;### at least 10 supporting reads

open(INDATA, $in_file)|| die "can not open file: \n";
<INDATA>; 
while($dataRaw=<INDATA>)
{
	$dataRaw=~s/\n|\r//g;
	if($dataRaw=~/\texon\t/)
	{
		@dataInf=split(/\s+/,$dataRaw);
		if($dataInf[4]>=$supporting_read)  ### at least 10 supporting reads
		{
			$dataInf[0]=~s/:|\|/_/g;
			$uniqueCircularRna{$dataInf[0]}++;
			@dataInf2=split(/\,/,$dataInf[9]);
			$geneWithCircularRna{$dataInf2[0]}++;
			$geneCircularRnaDetail{$dataInf2[0]}{$dataInf[0]}=0;
		}
	}		
}
close(INDATA);

print scalar(keys(%uniqueCircularRna)),"\n";
print scalar(keys(%geneWithCircularRna)),"\n";

open(IN, $ref_file);
open(OUT,">", $out_file);
while($raw=<IN>)
{
	if($raw=~/\texon\t/)
	{
		$raw1=$raw;
		$raw=~s/\n|\r|\"|\;//g;
		@inf=split(/\s+/,$raw);
		if(!exists($geneWithCircularRna{$inf[9]}))
		{
			print OUT $raw1;
		}
		else
		{
			print OUT $raw1;
			$hashMedian=$geneCircularRnaDetail{$inf[9]};
			%circularList=%$hashMedian;
			foreach $circular_location (keys(%circularList))
			{
				@circular_location_inf=split(/_/,$circular_location);
				if($inf[0] eq $circular_location_inf[0]&&$inf[3]>=$circular_location_inf[1]&&$inf[4]<=$circular_location_inf[2])
				{
					print OUT "$inf[0]\t$inf[1]\t$inf[2]\t$inf[3]\t$inf[4]\t$inf[5]\t$inf[6]\t$inf[7]\t$inf[8] \"$inf[9]\"\; $inf[10] \"circ_$circular_location\"\;\n";
				}
				if($inf[0] eq $circular_location_inf[0]&&$inf[3]<=$circular_location_inf[1]&&$inf[4]>=$circular_location_inf[2])
				{
					print OUT "$inf[0]\t$inf[1]\t$inf[2]\t$circular_location_inf[1]\t$circular_location_inf[2]\t$inf[5]\t$inf[6]\t$inf[7]\t$inf[8] \"$inf[9]\"\; $inf[10] \"circ_$circular_location\"\;\n";
				}
			}
		}
	}
	else
	{
		print OUT $raw;
	}
}

close(OUT);
close(IN);





















