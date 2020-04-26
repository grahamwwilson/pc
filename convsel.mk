# Makefile for the new ntuple processing
main: histset2.o convsel.o ParTreeProcessing2.C
	g++ -o compiledThreads2 ParTreeProcessing2.C -pthread `root-config --cflags --libs`

histset2.o: histset2.C convsel.o
	g++ -c -pthread histset2.C `root-config --cflags --libs`

convsel.o: convsel.C convsel.h
	g++ -c -pthread convsel.C `root-config --cflags --libs` 

clean:
	rm *.o
