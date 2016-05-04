#!/bin/bash

# Max for B is 2^7

echo "C|B|S|V|K|File Name|AAT"

for i in `ls ../traces/`; do
	C=15
	while [ $C -ge 0 ]; do
		if [ $C -gt 8 ]; then
			Bmax=8
		else
			Bmax=$C
		fi
		B=$Bmax
		while [ $B -ge 0 ]; do
			Smax=$((C - B))
			S=$Smax
			while [ $S -ge 0 ]; do
				V=4
				while [ $V -ge 0 ]; do
					K=4
					while [ $K -ge 0 ]; do
						./cachesim -c $C -b $B -s $S -v $V -k $K -i ../traces/$i		
						K=$((K - 1))  
					done;
					V=$((V - 1))  
				done;
				S=$((S - 1))  
			done;
			B=$((B - 1))  	  
		done;	
		C=$((C - 1))  
	done;
done;

#V=0
#while [ $V -lt 5 ]; do
#	K=0
#	while [ $K -lt 5 ]; do
#		for i in `ls ../traces/`; do
#				echo "--- $i ---"
#				echo ""
#				./cachesim -v $V -k $K -i ../traces/$i
#		done;
#		K=$((K + 1))
#	done;
#	V=$((V + 1))
#done;

