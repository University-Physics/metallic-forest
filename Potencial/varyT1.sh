arg=(1 2 5 7 10 15 20 30 40 50 70 100)
arg1=(1 2 3 4 5 6 7 8 9 10)
for T in ${arg[*]}; do
    for S in ${arg1[*]}; do
	echo $I
	echo $S
	./main.x $T 1 1 2 $S
    done
done
