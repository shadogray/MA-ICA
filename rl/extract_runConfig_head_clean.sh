#./extract_runConfig_head.sh |perl -pe 's/^(.*?;)(.*?;)/$2$1/' |sort -n |perl -p -e 's/;/ & /g; s/$/ \\\\/'
./extract_runConfig_head.sh |perl -pe 's/^(.*?;)(.*?;)/$2$1/' |sort -n |perl -p -e 's/;/ & /g; s/$/ \\\\/; s/^(\d.+?)(\d)/$1 \\textbf{$2}/; s/^(\d)/\\textbf{$1}/'
