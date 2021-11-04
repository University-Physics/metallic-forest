arg=(1 2 3 4 5 6 7 8 9 10)
arg1=(1 2 3 4 5 6 7 8 9 10)
arg2=(1 2 3 4 5 6 7 8 9 10)
for F in ${arg2[*]}; do
    for I in ${arg1[*]}; do
        for S in ${arg[*]}; do
	    ./main.x 1 1 1 0 $S $I $F &
	done
	wait
	for S in ${arg[*]}; do
	    mv "./data/Out1T1V1R0I"$I"S"$S"A"$F"F.txt" "./data/Frontera/Out1T1V1R0I"$I"S"$S"A"$F"F.txt"
	done
    done
done
