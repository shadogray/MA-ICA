cd $(dirname $0) || exit 1
for d in $*; do
	for i in 5; do 
		(cd ${d}/ && Rscript findWayFinLearn2.R . $i L) 
	done
done
