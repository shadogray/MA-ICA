pr_eval.sh ssteps[23456]/ |tee pr_eval.log 
head -1 correlation.csv.example > correlation.csv
cat pr_eval.log |./corr_extract.sh >> correlation.csv
Rscript correlation.R
