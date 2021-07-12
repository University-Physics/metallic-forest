arg=(6 7 8 9 10)
arg1=(9)
for I in ${arg1[*]}; do
    for S in ${arg[*]}; do
	echo $I
	echo $S
	./main.x 1  1 1 $I $S
	mv "./data/Out1T1V1R"$I"I"$S"S.txt" "./data/I1/Out1T1V1R"$I"I"$S"S.txt"
    done
done
