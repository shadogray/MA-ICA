cd $(dirname $0) || exit 1
DIRS=`echo $*|sed 's%/%%g'` 
for d in $DIRS; do
	for i in 2 3 4 5 6 7 8 9; do 
		(cd ${d}/ && Rscript findWayFinLearn2.R . $i LRp) 
	done
done
