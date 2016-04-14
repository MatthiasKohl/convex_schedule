CPP=g++
CFLAGS=-Ofast
LDFLAGS=
EXEC=primes
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

%.o: %.cpp
	$(CPP) -o $@ -c $< $(CFLAGS)

primes: primes.o
	$(CPP) $(CFLAGS) $(LDFLAGS) $(OBJ) -o $@

.PHONY: clean

clean:
	rm -rf *.o *.tmp $(EXEC)
