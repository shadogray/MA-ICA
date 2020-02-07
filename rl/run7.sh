cd $(dirname $0) || exit 1
for d in $*; do
	for i in 7; do 
		(cd ${d}/ && nohup nice -n 10 Rscript runFinLearn.R numSigs=${i} > runFinLearn_${i}.log &) 
	done
done
