CC = g++
CFLAGS = -Wall --std=c++11

VAR_GS.exe : VAR_GS.o funzioni.o random.o
	$(CC) random.o funzioni.o VAR_GS.o -o VAR_GS.exe
VAR_GS.o : VAR_GS.cpp VAR_GS.h
	$(CC) -c VAR_GS.cpp -o VAR_GS.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VAR_GS.exe seed.out
