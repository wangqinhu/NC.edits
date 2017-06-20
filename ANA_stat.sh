#!/bin/bash

#./extract_ANA.pl data/phy.red.txt > data/red.ANA.txt
#./extract_ANA.pl data/phy.rnd.txt > data/rnd.ANA.txt

ANA=("AAA" "AGA" "ACA" "ATA")

# nonsynonymous
echo "red_ANA_non"
for i in ${ANA[@]}; do
	grep non data/red.ANA.txt | cut -f3 | grep $i | wc -l | xargs echo $i
done
echo "rnd_ANA_non"
for i in ${ANA[@]}; do
	grep non data/rnd.ANA.txt | cut -f3 | grep $i | wc -l | xargs echo $i
done

# synonymous
echo "red_ANA_syn"
for i in ${ANA[@]}; do
	grep syn data/red.ANA.txt | cut -f3 | grep $i | wc -l | xargs echo $i
done
echo "rnd_ANA_syn"
for i in ${ANA[@]}; do
	grep syn data/rnd.ANA.txt | cut -f3 | grep $i | wc -l | xargs echo $i
done
