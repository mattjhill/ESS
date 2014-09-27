CC	=g++-4.9 -fopenmp
CFLAGS	=-c -Wall -Wno-write-strings -Ofast -ffast-math -ffinite-math-only -DNDEBUG -fomit-frame-pointer -I header/
LDFLAGS	=
SOURCES	=./testGRP.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./GRP

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out *.o GRP

tar:
	tar -zcvf GRP.tar.gz ./makefile.mk ./GRP.hpp ./testGRP.cpp ./README.md ./LICENSE.md
