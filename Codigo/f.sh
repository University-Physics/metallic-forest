arg1=(0 1 2 3)
for I in ${arg1[*]}; do
    ./aux.x 1 10 1 2 1 $I
done
