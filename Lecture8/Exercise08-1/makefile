CC = g++
CFLAGS = -g -c -Wall -std=c++11

SRCDIR = ../../random
OBJDIR = obj
RESDIR = results

EXECUTABLE = main
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

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