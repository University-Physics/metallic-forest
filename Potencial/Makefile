all:anim.gif
OBJ= funciones.o main.o
%.o:%.cpp constantes.h
	g++ -c $< -n

main.x:$(OBJ)

	g++ $^ -g -o $@

anim.gif: main.x
	./main.x|gnuplot
.PHONY:clean
clean:
	rm -f *.x *.o 
