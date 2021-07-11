arg=( 1 2 3 4 5 )
arg1=(1 2 3 4 5 6 7 8 9 10)
for R in ${arg[*]}; do
    for S in ${arg1[*]}; do
	echo $R
	echo $S
	./main.x 1 1 $R 2 $S
    done
done
