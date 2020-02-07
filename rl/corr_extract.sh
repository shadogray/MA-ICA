corr_eval.sh | sed -Ee 's%/% %' |perl -ne '/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*?S\.w\s+(.*)$/ && print "$1;$2;$3;$4;$5\n"' |sed 's/://'
