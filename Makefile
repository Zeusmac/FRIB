CC = g++
ROOTLIBS = `root-config --cflags --glibs`

all: Mainloop


Parcer: Parcer.cpp
	$(CC) -o Parcer.o -c Parcer.cpp $(ROOTLIBS)

TreeParcer: TreeParcer.cpp
	$(CC) -o TreeParcer TreeParcer.cpp $(ROOTLIBS)
	
Cloverdata: Cloverdata.cpp
	$(CC) -o Cloverdata Cloverdata.cpp $(ROOTLIBS)
	
YSO: YSO.cpp
	$(CC) -o YSO YSO.cpp $(ROOTLIBS)

OYSO: OYSO.cpp
	$(CC) -o OYSO.o -c YSO.cpp $(ROOTLIBS)

Clover: clover.cpp
	$(CC) -o clover.o -c clover.cpp $(ROOTLIBS)
	
Mainloop: Mainloop.cpp
	$(CC) -o Mainloop Mainloop.cpp $(ROOTLIBS)
	
Mainloopcatur: Mainloopcatur.cpp
	$(CC) -o Mainloopcatur Mainloopcatur.cpp $(ROOTLIBS)