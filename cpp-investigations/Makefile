CC := g++-8
STD := -std=c++17
CFLAGS := -O3 -fopenmp -mavx -march=native -funroll-loops

all: vectorize.cpp
	$(CC) $(STD) $(CFLAGS) vectorize.cpp -o vectorize
