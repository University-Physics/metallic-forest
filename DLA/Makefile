all:DLA.png
%.x:%.cpp
	g++ -std=c++11 $< -o $@
data.txt:DLA.x
	./$< > $@
DLA.png:grafica.py data.txt
	python3 $<
.PHONY:clean
clean:
	rm -f *.x
	rm DLA.png
