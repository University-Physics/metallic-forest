arg=(1 2 3 4 5 6 7 8 9 10)
arg1=(1 2 3 4 5 6 7 8 9 10)
for I in ${arg1[*]}; do
    for S in ${arg[*]}; do
	echo $I
	echo $S
	./main.x 100 100 1 $I $S
	mv "./data/Out100T100V1R"$I"I"$S"S.txt" "./data/I4/Out100T100V1R"$I"I"$S"S.txt"
    done
done
