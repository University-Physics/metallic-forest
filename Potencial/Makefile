all:anim.gif
OBJ= main.o
%.o:%.cpp constantes.h
	g++ -c -O2 $< -n

main.x:$(OBJ)

	g++ -O3 -g -Wall  $^  -o $@

anim.gif: main.x
	time ./main.x save movie | gnuplot
.PHONY:clean
clean:
	rm -f *.x *.o 
