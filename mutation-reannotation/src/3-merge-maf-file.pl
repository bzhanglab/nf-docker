#!usr/bin/perl -w

use Getopt::Long;

my ($input_gdc_maf, $input_maftools_maf, 
    $output_maf_full, $output_maf_2_tools, $output_maf_2_tools_or_hotspot,
    $mapping_file, $flag_genes_file);

Getopt::Long::GetOptions (
  'in1=s'          =>  \$input_gdc_maf,
  'in2=s'          =>  \$input_maftools_maf,
  'out1=s'          => \$output_maf_full,
  'out2=s'          => \$output_maf_2_tools,
  'out3=s'          => \$output_maf_2_tools_or_hotspot,
  'mapping_file=s'  => \$mapping_file,
  'flag_genes_file=s' => \$flag_genes_file
);

%flag_genes=();
open(IN,$flag_genes_file);
while($raw=<IN>)
{
    $raw=~s/\n|\r|\s+//g;
    $flag_genes{$raw}=1;
    #print "$raw\n";
}
close(IN);

%ensg_to_name=();
open(IN,$mapping_file);
while($raw=<IN>)
{
    $raw=~s/\n|\r//g;
    @inf=split(/\t/,$raw);
    $ensg_to_name{$inf[2]}=$inf[3];
}
close(IN);

%adopted_from_gdc=();
open(IN,$input_gdc_maf);
while($raw=<IN>)
{
    if($raw!~/^#/)
    {
        $raw=~s/\n|\r//g;
        @inf=split(/\t/,$raw);
        $adopted_from_gdc{"$inf[4]-$inf[5]"}="$inf[39]\t$inf[40]\t$inf[41]\t$inf[42]\t$inf[43]\t$inf[44]\t$inf[123]\t$inf[124]";
    }
}
close(IN);

open(IN,$input_maftools_maf);
<IN>;
open(OUT,">$output_maf_full");
print OUT "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification\ttx\texon\ttxChange\taaChange\tVariant_Type\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tavsnp150\tDamagePredCount\tSIFT_pred\tSIFT4G_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\tVEST4_score\tMetaSVM_pred\tMetaLR_pred\tM-CAP_pred\tREVEL_score\tMutPred_score\tMVP_score\tMPC_score\tPrimateAI_pred\tDEOGEN2_pred\tBayesDel_addAF_pred\tBayesDel_noAF_pred\tClinPred_pred\tLIST-S2_pred\tCADD_raw\tCADD_phred\tDANN_score\tfathmm-MKL_coding_pred\tfathmm-XF_coding_pred\tEigen-raw_coding\tEigen-phred_coding\tEigen-PC-raw_coding\tEigen-PC-phred_coding\tGenoCanyon_score\tintegrated_fitCons_score\tGM12878_fitCons_score\tH1-hESC_fitCons_score\tHUVEC_fitCons_score\tLINSIGHT\tGERP++_NR\tGERP++_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphyloP17way_primate\tphastCons100way_vertebrate\tphastCons30way_mammalian\tphastCons17way_primate\tbStatistic\tInterpro_domain\tGTEx_V8_gene\tGTEx_V8_tissue\tcosmic92Coding\tgene_name\tTumor_Sample_Barcode\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count\thotspot\tcallers\n";
while($raw=<IN>)
{
    $raw=~s/\n|\r//g;
    @inf=split(/\t/,$raw);
    if($raw=~/\tintergenic\t/)
    {
        $raw=~s/\n|\r//g;
        print OUT "$inf[6]\tBCM\tGRCh38.p13";
        for($i=1;$i<6;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        for($i=7;$i<=79;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        print OUT "\tNA";
        print OUT "\t$inf[0]";    
        print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";        
    }
    elsif(!exists($flag_genes{$ensg_to_name{$inf[6]}}))
    {
        $raw=~s/\n|\r//g;
        print OUT "$inf[6]\tBCM\tGRCh38.p13";
        for($i=1;$i<6;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        for($i=7;$i<=79;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        print OUT "\t$ensg_to_name{$inf[6]}";
        print OUT "\t$inf[0]";    
        print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";
    }
}
close(IN);
close(OUT);



open(IN,$input_maftools_maf);
<IN>;
open(OUT,">$output_maf_2_tools");
print OUT "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification\ttx\texon\ttxChange\taaChange\tVariant_Type\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tavsnp150\tDamagePredCount\tSIFT_pred\tSIFT4G_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\tVEST4_score\tMetaSVM_pred\tMetaLR_pred\tM-CAP_pred\tREVEL_score\tMutPred_score\tMVP_score\tMPC_score\tPrimateAI_pred\tDEOGEN2_pred\tBayesDel_addAF_pred\tBayesDel_noAF_pred\tClinPred_pred\tLIST-S2_pred\tCADD_raw\tCADD_phred\tDANN_score\tfathmm-MKL_coding_pred\tfathmm-XF_coding_pred\tEigen-raw_coding\tEigen-phred_coding\tEigen-PC-raw_coding\tEigen-PC-phred_coding\tGenoCanyon_score\tintegrated_fitCons_score\tGM12878_fitCons_score\tH1-hESC_fitCons_score\tHUVEC_fitCons_score\tLINSIGHT\tGERP++_NR\tGERP++_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphyloP17way_primate\tphastCons100way_vertebrate\tphastCons30way_mammalian\tphastCons17way_primate\tbStatistic\tInterpro_domain\tGTEx_V8_gene\tGTEx_V8_tissue\tcosmic92Coding\tgene_name\tTumor_Sample_Barcode\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count\thotspot\tcallers\n";
while($raw=<IN>)
{
    $raw=~s/\n|\r//g;
    @inf=split(/\t/,$raw);
    @inf_adopted=split("\t",$adopted_from_gdc{"$inf[1]-$inf[2]"});
    if($inf_adopted[@inf_adopted-1]=~/\;/)
    {
        if($raw=~/\tintergenic\t/)
        {
            $raw=~s/\n|\r//g;
            print OUT "$inf[6]\tBCM\tGRCh38.p13";
            for($i=1;$i<6;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            for($i=7;$i<=79;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            print OUT "\tNA";
            print OUT "\t$inf[0]";    
            print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";        
        }
        elsif(!exists($flag_genes{$ensg_to_name{$inf[6]}}))
        {
            $raw=~s/\n|\r//g;
            print OUT "$inf[6]\tBCM\tGRCh38.p13";
            for($i=1;$i<6;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            for($i=7;$i<=79;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            print OUT "\t$ensg_to_name{$inf[6]}";
            print OUT "\t$inf[0]";    
            print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";
        }
    }
}
close(IN);
close(OUT);



open(IN,$input_maftools_maf);
<IN>;
open(OUT,">$output_maf_2_tools_or_hotspot");
print OUT "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification\ttx\texon\ttxChange\taaChange\tVariant_Type\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tExAC_ALL\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tavsnp150\tDamagePredCount\tSIFT_pred\tSIFT4G_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\tVEST4_score\tMetaSVM_pred\tMetaLR_pred\tM-CAP_pred\tREVEL_score\tMutPred_score\tMVP_score\tMPC_score\tPrimateAI_pred\tDEOGEN2_pred\tBayesDel_addAF_pred\tBayesDel_noAF_pred\tClinPred_pred\tLIST-S2_pred\tCADD_raw\tCADD_phred\tDANN_score\tfathmm-MKL_coding_pred\tfathmm-XF_coding_pred\tEigen-raw_coding\tEigen-phred_coding\tEigen-PC-raw_coding\tEigen-PC-phred_coding\tGenoCanyon_score\tintegrated_fitCons_score\tGM12878_fitCons_score\tH1-hESC_fitCons_score\tHUVEC_fitCons_score\tLINSIGHT\tGERP++_NR\tGERP++_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphyloP17way_primate\tphastCons100way_vertebrate\tphastCons30way_mammalian\tphastCons17way_primate\tbStatistic\tInterpro_domain\tGTEx_V8_gene\tGTEx_V8_tissue\tcosmic92Coding\tgene_name\tTumor_Sample_Barcode\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count\thotspot\tcallers\n";
while($raw=<IN>)
{
    $raw=~s/\n|\r//g;
    @inf=split(/\t/,$raw);
    @inf_adopted=split("\t",$adopted_from_gdc{"$inf[1]-$inf[2]"});
    if($inf_adopted[@inf_adopted-2] eq "Y")
    {
        $raw=~s/\n|\r//g;
        print OUT "$inf[6]\tBCM\tGRCh38.p13";
        for($i=1;$i<6;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        for($i=7;$i<=79;$i++)
        {
            print OUT "\t$inf[$i]"
        }
        print OUT "\t$ensg_to_name{$inf[6]}";
        print OUT "\t$inf[0]";    
        print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";
    }    
    elsif($inf_adopted[@inf_adopted-1]=~/\;/)
    {
        if($raw=~/\tintergenic\t/)
        {
            $raw=~s/\n|\r//g;
            print OUT "$inf[6]\tBCM\tGRCh38.p13";
            for($i=1;$i<6;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            for($i=7;$i<=79;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            print OUT "\tNA";
            print OUT "\t$inf[0]";    
            print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";        
        }
        elsif(!exists($flag_genes{$ensg_to_name{$inf[6]}}))
        {
            $raw=~s/\n|\r//g;
            print OUT "$inf[6]\tBCM\tGRCh38.p13";
            for($i=1;$i<6;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            for($i=7;$i<=79;$i++)
            {
                print OUT "\t$inf[$i]"
            }
            print OUT "\t$ensg_to_name{$inf[6]}";
            print OUT "\t$inf[0]";    
            print OUT "\t",$adopted_from_gdc{"$inf[1]-$inf[2]"},"\n";
        }
    }
}
close(IN);
close(OUT);

