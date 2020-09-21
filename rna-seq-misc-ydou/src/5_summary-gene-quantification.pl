#!/usr/bin/perl -w

@sampleList=`ls * | grep isoforms.results`;

%geneRsem=();
%geneFpkm=();
%geneTpm=();
%geneIsoforms=();
%geneList=();
foreach $sampleName (@sampleList)
{
	$sampleName=~s/.isoforms.results//;
	print $sampleName;
	$sampleName=~s/\n|\r//g;
	open(INDATA,"$sampleName.genes.results")|| die "can not open file: $sampleName/$sampleName.genes.results\n";
	<INDATA>;
	while($dataRaw=<INDATA>)
	{
		@dataInf=();
		$dataRaw=~s/\n|\r//g;
		@dataInf=split(/\s+/,$dataRaw);
		$geneRsem{$sampleName}{$dataInf[0]}=$dataInf[4];
		$geneTpm{$sampleName}{$dataInf[0]}=$dataInf[5];
		$geneFpkm{$sampleName}{$dataInf[0]}=$dataInf[6];
		$geneIsoforms{$sampleName}{$dataInf[0]}=$dataInf[1];
		$geneList{$dataInf[0]}=1;
	}
	close(INDATA);
}

%circleIsoformRsem=();
%circleIsoformTpm=();
%circleIsoformFpkm=();
%circleIsoformRatio=();
%circleIsoformList=();
foreach $sampleName (@sampleList)
{
	$sampleName=~s/.isoforms.results//;
	print $sampleName,"\n";
	$sampleName=~s/\n|\r//g;
	open(INDATA,"$sampleName.isoforms.results")|| die "can not open file: $sampleName/$sampleName.isoforms.results\n";
	<INDATA>;
	while($dataRaw=<INDATA>)
	{
		@dataInf=();
		$dataRaw=~s/\n|\r//g;
		@dataInf=split(/\s+/,$dataRaw);
		if($dataInf[0]=~/circ_/)
		{
			$circleIsoformRsem{$sampleName}{"$dataInf[0]"."_"."$dataInf[1]"}=$dataInf[4];
			$circleIsoformTpm{$sampleName}{"$dataInf[0]"."_"."$dataInf[1]"}=$dataInf[5];
			$circleIsoformFpkm{$sampleName}{"$dataInf[0]"."_"."$dataInf[1]"}=$dataInf[6];
			$circleIsoformRatio{$sampleName}{"$dataInf[0]"."_"."$dataInf[1]"}=$dataInf[7];
			$circleIsoformList{"$dataInf[0]"."_"."$dataInf[1]"}=$dataInf[0];
		}
	}
	close(INDATA);
}
#################################################################################################################################
#################################start gene with circular rna#########################################################################
open(OUT,">gene_rsem_with_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneRsem{$sampleName}{$geneName}))
		{
			print OUT "\t$geneRsem{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">gene_fpkm_with_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneFpkm{$sampleName}{$geneName}))
		{
			print OUT "\t$geneFpkm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">gene_tpm_with_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneTpm{$sampleName}{$geneName}))
		{
			print OUT "\t$geneTpm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);
#################################################################################################################################
#################################end of gene with circular rna###################################################################
#################################################################################################################################
#################################start of gene without circular rna###################################################################

open(OUT,">gene_rsem_removed_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneRsem{$sampleName}{$geneName}))
		{
			$abundance=$geneRsem{$sampleName}{$geneName};
			@isoformList=split(/\,/,$geneIsoforms{$sampleName}{$geneName});
			foreach $isoformName (@isoformList)
			{
				if($isoformName=~/circ_/)
				{
					$abundance=$abundance-$circleIsoformRsem{$sampleName}{"$isoformName"."_"."$geneName"}
				}
			}
			print OUT "\t$abundance";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">gene_fpkm_removed_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneFpkm{$sampleName}{$geneName}))
		{
			$abundance=$geneFpkm{$sampleName}{$geneName};
			@isoformList=split(/\,/,$geneIsoforms{$sampleName}{$geneName});
			foreach $isoformName (@isoformList)
			{
				if($isoformName=~/circ_/)
				{
					$abundance=$abundance-$circleIsoformFpkm{$sampleName}{"$isoformName"."_"."$geneName"}
				}
			}
			print OUT "\t$abundance";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">gene_tpm_removed_circRNA_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneTpm{$sampleName}{$geneName}))
		{
			$abundance=$geneTpm{$sampleName}{$geneName};
			@isoformList=split(/\,/,$geneIsoforms{$sampleName}{$geneName});
			foreach $isoformName (@isoformList)
			{
				if($isoformName=~/circ_/)
				{
					$abundance=$abundance-$circleIsoformTpm{$sampleName}{"$isoformName"."_"."$geneName"}
				}
			}
			print OUT "\t$abundance";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

#################################################################################################################################
#################################end of gene without circular rna###################################################################
#################################################################################################################################
#################################start of circular rna###################################################################

open(OUT,">circRNA_rsem_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%circleIsoformList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($circleIsoformRsem{$sampleName}{$geneName}))
		{
			print OUT "\t$circleIsoformRsem{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">circRNA_fpkm_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%circleIsoformList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($circleIsoformFpkm{$sampleName}{$geneName}))
		{
			print OUT "\t$circleIsoformFpkm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);


open(OUT,">circRNA_tpm_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%circleIsoformList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($circleIsoformTpm{$sampleName}{$geneName}))
		{
			print OUT "\t$circleIsoformTpm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

#################################################################################################################################
#################################isoform level ###################################################################
#################################################################################################################################
@sampleList=`ls * | grep isoforms.results`;

%geneRsem=();
%geneFpkm=();
%geneTpm=();
%geneIsoforms=();
%geneList=();
foreach $sampleName (@sampleList)
{
	$sampleName=~s/.isoforms.results//;
	print $sampleName;
	$sampleName=~s/\n|\r//g;
	open(INDATA,"$sampleName.isoforms.results")|| die "$sampleName.isoforms.results\n";
	<INDATA>;
	while($dataRaw=<INDATA>)
	{
		if($dataRaw!~/circ_/)
		{
			@dataInf=();
			$dataRaw=~s/\n|\r//g;
			@dataInf=split(/\s+/,$dataRaw);
			$geneRsem{$sampleName}{$dataInf[0]}=$dataInf[4];
			$geneTpm{$sampleName}{$dataInf[0]}=$dataInf[5];
			$geneFpkm{$sampleName}{$dataInf[0]}=$dataInf[6];
			$geneIsoforms{$sampleName}{$dataInf[0]}=$dataInf[1];
			$geneList{$dataInf[0]}=1;
		}
	}
	close(INDATA);
}

#################################################################################################################################
open(OUT,">linear_isoform_rsem_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneRsem{$sampleName}{$geneName}))
		{
			print OUT "\t$geneRsem{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">linear_isoform_fpkm_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneFpkm{$sampleName}{$geneName}))
		{
			print OUT "\t$geneFpkm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">linear_isoform_tpm_BCM.txt");
print OUT "idx";
foreach $sampleName (@sampleList)
{
	print OUT "\t$sampleName";
}
print OUT "\n";
foreach $geneName (sort keys(%geneList))
{
	print OUT "$geneName";
	foreach $sampleName (@sampleList)
	{
		if(exists($geneTpm{$sampleName}{$geneName}))
		{
			print OUT "\t$geneTpm{$sampleName}{$geneName}";
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);
