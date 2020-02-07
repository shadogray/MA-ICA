#!/bin/bash
egrep ' rowname|Correlation|^. S.w' | perl -0777 -ne 'while(/(ICA.*?) (.*?)\n(.*?)\n(\d S\.w)\s+(-?)(.*?) (.*?)\n/gs) { print "$6 $1 $5$6 $7 $1 $4  $2\n"}'
