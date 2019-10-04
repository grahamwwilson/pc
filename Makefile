main: histset.o myselector.o ParTreeProcessing.C
	g++ -o compiledthreads ParTreeProcessing.C -pthread `root-config --cflags --libs`

histset.o: histset.C myselector.o
	g++ -c  histset.C -pthread `root-config --cflags --libs`

myselector.o: myselector.C myselector.h
	g++ -c -pthread myselector.C `root-config --cflags --libs` 

clean:
	rm *.o
