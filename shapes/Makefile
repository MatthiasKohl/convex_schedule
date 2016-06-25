CPP=g++
CFLAGS=-Ofast
LDFLAGS=
EXEC=gen_shapes
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

%.o: %.cpp
	$(CPP) -o $@ -c $< $(CFLAGS)

gen_shapes: gen_shapes.o
	$(CPP) $(CFLAGS) $(LDFLAGS) $(OBJ) -o $@

.PHONY: clean

clean:
	rm -rf *.o *.tmp $(EXEC)
