find . -type f -empty -print -delete
for f in *.score.txt; 
do 
LINECOUNT=`wc -l $f | cut -f1 -d ' '`; 


if [[ $LINECOUNT -eq 1 ]]; then
rm $f
echo $f $LINECOUNT; 
fi

done


for f in ~/scratch/ENCO*_*/output/variants/clinvar;
do
cd $f
echo $f
pwd
find . -name '*.vep.tsv' -type f -empty -print 
done