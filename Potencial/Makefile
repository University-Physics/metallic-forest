all:main.x
OBJ= main.o

%.o:%.cpp constantes.h
	g++ -c -O3 $< -n

main.x:$(OBJ)
	g++ -O3 -g -pg -Wall $^ -o $@

I1: main.x varyI1.sh
	nohup bash varyI1.sh

I2: main.x varyI2.sh
	nohup bash varyI2.sh

I3: main.x varyI3.sh
	nohup bash varyI3.sh

I4: main.x varyI4.sh
	nohup bash varyI4.sh

R1: main.x varyR1.sh
	nohup bash varyR1.sh

R2: main.x varyR2.sh
	nohup bash varyR2.sh

R3: main.x varyR3.sh
	nohup bash varyR3.sh

R4: main.x varyR4.sh
	nohup bash varyR4.sh

T1: main.x varyT1.sh
	nohup bash varyT1.sh

T2: main.x varyT2.sh
	nohup bash varyI2.sh

V1: main.x varyV1.sh
	nohup bash varyV1.sh

V2: main.x varyV2.sh
	nohup bash varyI1.sh

Fractal_imbalance.png: main.x run.sh Fractal_dim.py
	bash run.sh
	python3 Fractal_dim.py

.PHONY:clean
clean:
	rm -f *.x *.o 
