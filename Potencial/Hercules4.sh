arg=(11 12 13 14 15 16 17 18 19 20)
arg1=(1 2 3 4 5 )
for V in ${arg[*]}; do
    for S in ${arg1[*]}; do
	echo $I
	echo $S
	./main.x 1 $V 1 2 $S &
    done
    wait
    for S in ${arg1[*]}; do
	mv "./data/Out1T"$V"V1R2I"$S"S.txt" "./data/V_examine/Out1T"$V"V1R2I"$S"S.txt"
    done
done
