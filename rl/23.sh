cd $(dirname $0) || exit 1
for i in 2 3; do (cd $1/ && nohup nice -n 10 Rscript runFinLearn${i}.R > runFinLearn_ica_${i}.log &) done
