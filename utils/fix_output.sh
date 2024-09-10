find . -type f -empty -print -delete
for f in *.score.txt; 
do 
LINECOUNT=`wc -l $f | cut -f1 -d ' '`; 


if [[ $LINECOUNT -eq 1 ]]; then
rm $f
echo $f $LINECOUNT; 
fi

done