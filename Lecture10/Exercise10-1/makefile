CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o punto.o
	$(CC) random.o punto.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
punto.o : punto.cpp punto.h
	$(CC) -c punto.cpp -o punto.o $(CFLAGS)
random.o : ../../random/random.cpp ../../random/random.h
	$(CC) -c ../../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove :
	rm *.txt