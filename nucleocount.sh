#!/bin/bash
TEST=$(cat "$1")
for i in $TEST; do
	java -Xmx1024m  -Djava.awt.headless=true -jar ~/IGVTools/igvtools.jar count -w 1 --bases --query $i $2 count.wig $3
	BASES=$(sed -n '4p' count.wig)
	OUTTEXT=$i$'\t'$BASES
	echo $OUTTEXT >> $4
done
