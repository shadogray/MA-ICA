FILES=`ls ssteps[0-9]/*[0-9].rdata |sed -e 's/finDemo.env_monster_ica_/runFinLearn/' -e 's/rdata/R/'` 

for f in $FILES; do grep -v deltaSumR $f | perl -0777 -ple 's/\n/;/gm; s/ +/ /g; s/^.*?(numSig.*?gamma *= *[\d.]+).*$/\1/; s/;+/;/g'; echo ""; done |perl -pe 's/(\w+) = [*\d.]+/\1/g' |head -1

for f in $FILES; do grep -v deltaSumR $f | perl -0777 -ple 's/\n/;/gm; s/ +/ /g; s/^.*?(numSig.*?gamma *= *[\d.]+).*$/\1/; s/;+/;/g'; echo ""; done |perl -pe 's/\w+ = ([*\d.]+)/\1/g' |sort -n

