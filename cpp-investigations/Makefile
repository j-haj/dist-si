CC := g++-7
STD := -std=c++17
CFLAGS := -Ofast -fopenmp -mavx -march=native -funroll-loops

all: vectorize.cpp
	$(CC) $(STD) $(CFLAGS) vectorize.cpp -o vectorize
