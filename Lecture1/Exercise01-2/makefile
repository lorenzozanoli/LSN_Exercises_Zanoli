CC = g++
CFLAGS = -g -c -Wall -std=c++11

SRCDIR = ../../random
OBJDIR = obj
RESDIR = results

EXECUTABLE = main
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

exe: $(EXECUTABLE)
	./main

$(EXECUTABLE): $(OBJECTS) $(OBJDIR)/main.o
	$(CC) $(LFLAGS) $^ -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) -o $@ $<

$(OBJDIR)/main.o: main.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) -o $@ $<

$(OBJDIR):
	mkdir $(OBJDIR)

clean:
	rm -f -r $(OBJDIR) *.o $(EXECUTABLE)

cleanall:
	rm -f -r $(OBJDIR) $(RESDIR) *.o $(EXECUTABLE)


# ISTRUZIONI PER L'USO \
-> La funzione $(patsubst PATTERN,REPLACEMENT,TEXT) prende le parole in TEXT \
   compatibili con PATTERN e le sostituisce con REPLACEMENT. Nel mio caso prende \
   in RANNYU tutti i nomi di file .cc e li sostituisce con .o \
   Non fa nulla di compilazione, serve solo per scrivere i target negli argomenti. \
-> $@ riprende il nome della regola \
   $< riprende il nome del primo argomento \
   $^ riprende tutti i prerequisiti