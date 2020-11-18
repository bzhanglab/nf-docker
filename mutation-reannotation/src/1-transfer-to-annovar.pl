#!usr/bin/perl -w
use Getopt::Long;

my ($sample_id, $in, $out, $help);

Getopt::Long::GetOptions (
  'in|I=s'          =>  \$in,
  'out|O=s'         =>  \$out,
  'help|H!'         =>  \$help
);


if (defined($help)) {
  print "
Usage: ./1-transfer-to-annovar.pl -I in.maf -O output.txt

Arguments:
    -I, --in
          input maf file
    -O, --out
          output converted file
\n";
    exit 0;
}

open(IN, $in);
open(OUT, ">", $out);
while($raw=<IN>)
{
    if($raw!~/^#/&&$raw!~/^Hugo_Symbol/)
    {
        $raw=~s/\n|\r//g;
        @inf=split(/\t/,$raw);
        if($inf[11] eq "-") ###deletions
        {
            print OUT "$inf[4]\t$inf[5]\t$inf[5]\t$inf[11]\t$inf[12]\n";
        }
        else
        {
            print OUT "$inf[4]\t$inf[5]\t$inf[6]\t$inf[11]\t$inf[12]\n";
        }
    }
}
close(IN);
close(OUT);
