all:anim.gif
OBJ= main.o
%.o:%.cpp constantes.h
	g++ -c -O3 -fopenmp $< -n

main.x:$(OBJ)

	g++ -O3 -g -pg -Wall -fopenmp $^  -o $@

anim.gif: main.x
	time ./main.x save movie | gnuplot
	python3 grafica.py
.PHONY:clean
clean:
	rm -f *.x *.o 
