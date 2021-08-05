arg=(1 2 3 4 5 6 7 8 9 10)
for V in ${arg[*]}; do
    ./main.x 1 10 1 2 $V $1 &
done
wait
for S in ${arg[*]}; do
    mv "./data/Out1T10V1R2I"$S"S.txt" "./data/fronteras/"$1"Out1T10V1R2I"$S"S.txt"
    mv "1T10V1R2I"$S"SProbability_distribution.txt" "./data/fronteras/P/1T10V1R2I"$S"SProbability_distribution.txt"
done

