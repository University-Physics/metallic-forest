arg1=(1 2 3 4 5 6 7 8 9 10)
for S in ${arg1[*]}; do
    ./main.x 100 1 1 2 $S 1
done
