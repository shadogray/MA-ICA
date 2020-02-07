cd $(dirname $0) || exit 1
for d in $*; do
	(for i in 2 3 4 5; do 
		(cd ${d}/ && Rscript findWayFinLearn2.R . $i LpR) 
	done) |tee print_${d}.log 
done
