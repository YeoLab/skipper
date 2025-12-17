mkdir ~/ps-yeolab5/roulette
cd ~/ps-yeolab5/roulette
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
wget http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${c}_rate_v5.2_TFBS_correction_all.vcf.bgz
wget http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${c}_rate_v5.2_TFBS_correction_all.vcf.bgz.csi
done