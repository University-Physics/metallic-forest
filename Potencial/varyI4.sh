arg=(1 2 3 4 5 6 7 8 9 10)
arg1=(1 2 3 4 5 6 7 8 9 10)
for I in ${arg1[*]}; do
    for S in ${arg[*]}; do
	echo $I
	echo $S
	./main.x 100 100 1 $I $S
    done
done
