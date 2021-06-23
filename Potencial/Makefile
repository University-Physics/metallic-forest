all:Fractal_imbalance.png
OBJ= main.o

%.o:%.cpp constantes.h
	g++ -c -O3 $< -n

main.x:$(OBJ)
	g++ -O3 -g -pg -Wall $^ -o $@

Fractal_imbalance.png: main.x run.sh Fractal_dim.py
	bash run.sh
	python3 Fractal_dim.py

.PHONY:clean
clean:
	rm -f *.x *.o 
