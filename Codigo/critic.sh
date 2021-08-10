arg1=(0 1 2 3 4 5 6 7 8 9 10)
for I in ${arg1[*]}; do
	./main.x 1 2 1 0 $I
done

