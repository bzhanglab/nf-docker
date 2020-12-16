#!/bin/bash

declare -A ASSEMBLIES
ASSEMBLIES[hs37d5]="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
ASSEMBLIES[hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz"
ASSEMBLIES[GRCh37]="ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
ASSEMBLIES[hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromFa.tar.gz"
ASSEMBLIES[GRCh38]="ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
ASSEMBLIES[mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz"
ASSEMBLIES[GRCm38]="ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
ASSEMBLIES[GRCh38.p13]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz"

declare -A ANNOTATIONS
ANNOTATIONS[GENCODE19]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL87]="ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
ANNOTATIONS[GENCODE28]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL93]="ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.chr.gtf.gz"
ANNOTATIONS[GENCODEM25]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
ANNOTATIONS[RefSeq_mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz"
ANNOTATIONS[GENCODEV34_basic]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.basic.annotation.gtf.gz"

declare -A COMBINATIONS
COMBINATIONS["hs37d5+GENCODE19"]="hs37d5+GENCODE19"
COMBINATIONS["hs37d5+RefSeq"]="hs37d5+RefSeq_hg19"
COMBINATIONS["hs37d5+ENSEMBL87"]="hs37d5+ENSEMBL87"
COMBINATIONS["hg19+GENCODE19"]="hg19+GENCODE19"
COMBINATIONS["hg19+RefSeq"]="hg19+RefSeq_hg19"
COMBINATIONS["hg19+ENSEMBL87"]="hg19+ENSEMBL87"
COMBINATIONS["GRCh37+GENCODE19"]="GRCh37+GENCODE19"
COMBINATIONS["GRCh37+RefSeq"]="GRCh37+RefSeq_hg19"
COMBINATIONS["GRCh37+ENSEMBL87"]="GRCh37+ENSEMBL87"
COMBINATIONS["hg38+GENCODE28"]="hg38+GENCODE28"
COMBINATIONS["hg38+RefSeq"]="hg38+RefSeq_hg38"
COMBINATIONS["hg38+ENSEMBL93"]="hg38+ENSEMBL93"
COMBINATIONS["GRCh38+GENCODE28"]="GRCh38+GENCODE28"
COMBINATIONS["GRCh38+RefSeq"]="GRCh38+RefSeq_hg38"
COMBINATIONS["GRCh38+ENSEMBL93"]="GRCh38+ENSEMBL93"
COMBINATIONS["GRCm38+GENCODEM25"]="GRCm38+GENCODEM25"
COMBINATIONS["GRCh38.p13+GENCODEV34_basic"]="GRCh38.p13+GENCODEV34_basic"
COMBINATIONS["GRCm38+RefSeq"]="GRCm38+RefSeq_mm10"
COMBINATIONS["mm10+GENCODEM25"]="mm10+GENCODEM25"
COMBINATIONS["mm10+RefSeq"]="mm10+RefSeq_mm10"
for COMBINATION in ${!COMBINATIONS[@]}; do
	COMBINATIONS["${COMBINATION%+*}viral+${COMBINATION#*+}"]="${COMBINATIONS[$COMBINATION]%+*}viral+${COMBINATIONS[$COMBINATION]#*+}"
done

if [ $# -ne 1 ] || [ -z "$1" ] || [ -z "${COMBINATIONS[$1]}" ]; then
	echo "Usage: $(basename $0) ASSEMBLY+ANNOTATION" 1>&2
	echo "Available assemblies and annotations:" 1>&2
	sed -e 's/ /\n/g' <<<"${!COMBINATIONS[@]}" | sort 1>&2
	exit 1
fi

ASSEMBLY="${COMBINATIONS[$1]%+*}"
VIRAL=""
if [[ $ASSEMBLY =~ viral ]]; then
	ASSEMBLY="${ASSEMBLY%viral}"
	VIRAL="viral"
fi
ANNOTATION="${COMBINATIONS[$1]#*+}"
THREADS="${THREADS-8}"
SJDBOVERHANG="${SJDBOVERHANG-250}"

set -o pipefail
set -e -u

echo "Downloading assembly: ${ASSEMBLIES[$ASSEMBLY]}"
wget -q -O - "${ASSEMBLIES[$ASSEMBLY]}" |
if [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.tar\.gz$ ]]; then
	tar -x -O -z
elif [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.gz$ ]]; then
	gunzip -c
else
	cat
fi |
if [ "$VIRAL" = "viral" ]; then
	# drop viral contigs from hs37d5 assembly
	awk '/^>/{ contig=$1 } contig!~/^>NC_|^>AC_/{ print }'
else
	cat
fi > "$ASSEMBLY$VIRAL.fa"

if [ "$VIRAL" = "viral" ]; then
	echo "Appending RefSeq viral genomes"
	REFSEQ_VIRAL_GENOMES=$(dirname "$0")/RefSeq_viral_genomes_v2.0.0.fa.gz
	if [ ! -e "$REFSEQ_VIRAL_GENOMES" ]; then
		REFSEQ_VIRAL_GENOMES=$(dirname "$0")/database/RefSeq_viral_genomes_v2.0.0.fa.gz
	fi
	gunzip -c "$REFSEQ_VIRAL_GENOMES" >> "$ASSEMBLY$VIRAL.fa"
fi

echo "Downloading annotation: ${ANNOTATIONS[$ANNOTATION]}"
wget -q -O - "${ANNOTATIONS[$ANNOTATION]}" |
if [[ ${ANNOTATIONS[$ANNOTATION]} =~ \.gz$ ]]; then
	gunzip -c
else
	cat
fi |
if [[ $ANNOTATION =~ RefSeq ]]; then
	# convert genePred to GTF
	awk -F '\t' -v OFS='\t' '
	function min(x, y) { return (x>y) ? y : x }
	function max(x, y) { return (x<y) ? y : x }
	{
		split($10, start, ",")
		split($11, end, ",")
		split($16, frame, ",")
		# remove stop codon from left end for coding genes on the minus strand
		if ($4=="-" && $14=="cmpl" && (start[1]!=$7 || (min(end[1],$8)-start[1]+frame[1])%3==0)) {
			$7+=3
			for (i in end)
				if ($7>=end[i] && $7<=end[i]+2)
					$7+=start[i+1]-end[i]
		}
		# remove stop codon from right end for coding genes on the plus strand
		if ($4=="+" && $15=="cmpl" && (end[$9]!=$8 || (end[$9]-max(start[$9],$7)+frame[$9])%3==0)) {
			$8-=3
			for (i in start)
				if ($8<=start[i] && $8>=start[i]-2)
					$8-=start[i]-end[i-1]
		}
		# append running number to duplicate uses of the same transcript ID
		gene_id=$13
		if (transcripts[$2]++) {
			gene_id=$13"_"transcripts[$2]
			$2=$2"_"transcripts[$2]
		}
		# print one line for each exon
		for (i=1; i<=$9; i++) {
			exon=($4=="+") ? i : $9-i+1
			attributes="gene_id \""gene_id"\"; transcript_id \""$2"\"; exon_number \""exon"\"; exon_id \""$2"."exon"\"; gene_name \""$13"\";"
			print $3,"RefSeq","exon",start[i]+1,end[i],".",$4,".",attributes
			# print one line for each coding region
			if ($14~/cmpl/ && $7<=end[i] && $8>=start[i])
				print $3,"RefSeq","CDS",max($7,start[i])+1,min($8,end[i]),".",$4,frame[i],attributes
		}
	}' | sort -k1,1V -k4,4n -k5,5n -k3,3 -S4G
else
	cat
fi |
if ! grep -q '^>chr' "$ASSEMBLY$VIRAL.fa"; then
	sed -e 's/^chrM\t/MT\t/' -e 's/^chr//'
else
	sed -e 's/^MT\t/chrM\t/' -e 's/^\([1-9XY]\|[12][0-9]\)\t/\1\t/'
fi > "$ANNOTATION.gtf"

mkdir STAR_index_${ASSEMBLY}${VIRAL}_${ANNOTATION}
STAR --runMode genomeGenerate --genomeDir STAR_index_${ASSEMBLY}${VIRAL}_${ANNOTATION} --genomeFastaFiles "$ASSEMBLY$VIRAL.fa" --sjdbGTFfile "$ANNOTATION.gtf" --runThreadN "$THREADS" --sjdbOverhang "$SJDBOVERHANG"
tar cvfz STAR_index_${ASSEMBLY}${VIRAL}_${ANNOTATION}.tgz STAR_index_${ASSEMBLY}${VIRAL}_${ANNOTATION}

