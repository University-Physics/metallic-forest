arg=(1 2 3 4 5 6 7 8 9 10)
arg1=(1 2 3 4 5 6 7 8 9 10)
for V in ${arg[*]}; do
    for S in ${arg1[*]}; do
	mv "./data/V2/Out100T"$V"V1R2I"$S"S.txt" "./data/V2/1Out100T"$V"V1R2I"$S"S.txt"
    done
done