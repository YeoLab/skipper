# non-RBP site introns
intron=/tscc/nfs/home/hsher/gencode_coords/gencode.v35.basic.intron.gff3

# bedtools intersect multiple rbps # rename stupid shit
cat ~/scratch/ENCODE3_HepG2/output/finemapping/mapped_sites/*finemapped_windows.bed.gz > ~/scratch/allrbpsite.bed.gz

bedtools intersect -a $intron -b ~/scratch/allrbpsite.bed.gz \
-s -v | awk '{ gsub(/chr/,"", $1); print }' OFS="\t" | \
awk '{ print $1,$4,$5,$6,$3,$7}'  OFS="\t" > ~/scratch/non_rbp_intron.bed

# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

bcftools query -R ~/scratch/non_rbp_intron.bed \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
    /tscc/nfs/home/hsher/ps-yeolab5/roulette/22_rate_v5.2_TFBS_correction_all.vcf.bgz \
    | awk '$1="chr"$1' OFS="\t" > ~/scratch/non_rbp_intron.roulette.vcf

# rename
awk '$1="chr"$1' OFS="\t" ~/scratch/non_rbp_intron.bed > ~/scratch/non_rbp_intron.chr.bed


# query genes
gencode=/tscc/nfs/home/hsher/gencode_coords/gencode.v40.primary_assembly.annotation.gff3
genes_bed=~/scratch/genes.bed
gene_vcf=~/scratch/genes.vcf

awk '$3=="gene" {print $1,$4,$5,$6,$8,$7}' OFS="\t" $gencode > $genes_bed

bcftools query -R $genes_bed \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
    /tscc/nfs/home/hsher/ps-yeolab5/gnomAD/v3/gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz \
    > $gene_vcf


bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
    /tscc/nfs/home/hsher/ps-yeolab5/gnomAD/v3/gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz \
    > ~/scratch/gnomad.chr22.tsv

bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/MR\t%INFO/AR\t%INFO/MG\t%INFO/MC\n' \
    /tscc/nfs/home/hsher/ps-yeolab5/roulette/22_rate_v5.2_TFBS_correction_all.header.filtered.vcf \
    > ~/scratch/roulette.chr22.vcf
