CC = g++
CFLAGS = -O3 -std=c++11

NVE.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE.exe

MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)

random.o : ../../random/random.cpp ../../random/random.h
	$(CC) -c ../../random/random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.exe

cleanall:
	rm *.o *.exe
	rm /output/solid/*.dat
	rm /output/liquid/*.dat
	rm /output/gas/*.dat