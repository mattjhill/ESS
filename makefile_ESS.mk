CC		=g++
CFLAGS	=-c -Wall -Ofast -ffast-math -ffinite-math-only -DNDEBUG -I header/
LDFLAGS	=
SOURCES	=./testESS.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./ESS

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out *.o ESS

tar:
	tar -zcvf ESS.tar.gz ./makefile.mk ./ESS.hpp ./testESS.cpp ./README.md ./LICENSE.md
