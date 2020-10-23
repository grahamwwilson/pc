# Makefile for the new ntuple processing
main: histset3.o recosim.o ParTreeProcessing3.C
	g++ -o compiledThreads3 ParTreeProcessing3.C -pthread `root-config --cflags --libs`

histset3.o: histset3.C recosim.o Hungarian.h
	g++ -c -pthread histset3.C `root-config --cflags --libs`

recosim.o: recosim.C recosim.h
	g++ -c -pthread recosim.C `root-config --cflags --libs`

hung.o: Hungarian.cpp Hungarian.h
	g++ -c -pthread Hungarian.cpp `root-config --cflags --libs`

clean:
	rm *.o
