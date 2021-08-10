arg=(1 2 5 7 10 15 20 30 40 50 70 100 200 250 300 350 400 500 600 700)
arg1=(1 2 3 4 5 6 7 8 9 10)
for V in ${arg[*]}; do
    for S in ${arg1[*]}; do
	echo $I
	echo $S
	./main.x 1 $V 1 2 $S
	mv "./data/Out1T"$V"V1R2I"$S"S.txt" "./data/V1/Out1T"$V"V1R2I"$S"S.txt"
    done
done
