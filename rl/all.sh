cd $(dirname $0) || exit 1
for d in $*; do
	for i in 2 3 4 5 6 7 8; do 
		(cd ${d}/ && nohup nice -n 10 Rscript runFinLearn${i}.R > runFinLearn_ica_${i}.log &) 
	done
done
