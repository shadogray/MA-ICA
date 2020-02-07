cd $(dirname $0) || exit 1
DIRS=`echo $*|sed 's%/%%g'` 
for d in $DIRS; do
	(
	for i in 2 3 4 5 6 7 8 9; do 
		#(cd ${d}/ && xvfb-run Rscript findWayFinLearn2.R . $i L) 
		(cd ${d}/ && Rscript findWayFinLearn2.R . $i L) 
	done
	tar cfv png_${d}.tar ${d}/*.png 
	tar cfv rdata_${d}.tar ${d}/*_Res.rdata
	echo "done: $d"
	)&
done
