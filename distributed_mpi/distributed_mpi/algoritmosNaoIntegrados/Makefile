MAIN = mpi_medians.c

TARGET = ballAlg

FLAGS = -O3 -fopenmp -lm
debugFLAGS = -fopenmp -g -lm -Wall

all: ballAlg

$(TARGET): $(MAIN) gen_points.o psrs_mpi.o
	mpicc $(debugFLAGS)  $^ -o $@ $(debugFLAGS)

gen_points.o: gen_points.c
	mpicc $^ -c $(debugFLAGS)

psrs_mpi.o: psrs_mpi.c
	mpicc $^ -c $(debugFLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
