MAIN = ballAlg-omp.c

TARGET = ballAlg-omp

OBJS = 

FLAGS = -O3 -fopenmp -lm
debugFLAGS = -fopenmp -g -lm -Wall

all: ballAlg-omp

$(TARGET): gen_points.o geometry.o unsorted_median.o
	gcc $(MAIN)  $^ -o $@ $(FLAGS)

gen_points.o: gen_points.c
	gcc $^ -c $(FLAGS)

geometry.o: geometry.c
	gcc $^ -c $(FLAGS)

unsorted_median.o:unsorted_median.c
	gcc $^ -c $(FLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
