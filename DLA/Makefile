all:DLA.png
OBJ= functions-DLA.o main-DLA.o
%.o:%.cpp
	g++ -c $< -n

main-DLA.x:$(OBJ)	
	g++ $^ -o $@

data.txt:main-DLA.x 
	./$< > $@
DLA.png:grafica.py data.txt
	python3 $<
.PHONY:clean
clean:
	rm -f *.x
	rm -f *.o
	rm -f *.txt
	rm DLA.png
