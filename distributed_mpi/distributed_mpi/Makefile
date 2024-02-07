MAIN = ballAlg-mpi.c

TARGET = ballAlg

OBJS = 

FLAGS = -O3 -fopenmp -lm
debugFLAGS = -fopenmp -g -lm -Wall

all: ballAlg

$(TARGET): gen_points.o geometry.o unsorted_median.o
	mpicc $(MAIN)  $^ -o $@ $(FLAGS)

gen_points.o: gen_points.c
	mpicc $^ -c $(FLAGS)

geometry.o: geometry.c
	mpicc $^ -c $(FLAGS)

unsorted_median.o:unsorted_median.c
	mpicc $^ -c $(FLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
