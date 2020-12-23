all :  lineareq

lineareq : main.o input.o inverse.o 
	mpicc main.o input.o inverse.o -lm -o ./a.out   -O3

main.o : main.c inverse.h input.h
	mpicc -c main.c   		-O3
	
input.o : input.c input.h
	mpicc -c input.c 		-O3

inverse.o : inverse.c  inverse.h
	mpicc -c inverse.c 		-O3
	
clean:
	rm -rf *.o ./a.out

#.o - объектные файлы
# Запуск : mpirun  -np <число запусков программы> ./a.out ...
