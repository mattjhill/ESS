CC		=g++
CFLAGS	=-c -Wall -Wno-write-strings -Ofast -ffast-math -ffinite-math-only -DNDEBUG -fomit-frame-pointer -I header/
LDFLAGS	=
SOURCES	=./testGRPpaperResults.cpp
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./GRPpaperResults

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.out *.o GRPpaperResults

tar:
	tar -zcvf GRP.tar.gz ./makefile.mk ./GRP.hpp ./testGRP.cpp ./testGRPpaperResults.cpp ./README.md ./LICENSE.md
