all:main.x
OBJ= main.o
i=1
r=1
t=1
v=1

%.o:%.cpp constantes.h
	g++ -c -O3 $< 

main.x:$(OBJ)
	g++ -O3 -g -Wall -fsanitize=undefined -fsanitize=address $^ -o $@

I$(i): main.x varyI$(i).sh
	nohup bash varyI$(i).sh

R$(r): main.x varyR$(r).sh
	nohup bash varyR$(r).sh

T$(t): main.x varyT$(t).sh
	nohup bash varyT$(t).sh

V$(v): main.x varyV$(v).sh
	nohup bash varyV$(v).sh

Fractal_imbalance.png: main.x run.sh Fractal_dim.py
	bash run.sh
	python3 Fractal_dim.py
aux.x:aux.cpp
	g++ -O3 -g -Wall $< -o $@
.PHONY:clean
clean:
	rm -f *.x *.o 
