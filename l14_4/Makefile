
prog_name = l14_4

CC = g++ -g -std=c++11
#CC = gcc -g



all: clean comp run

clean:
	rm -rf ./${prog_name}.exe

comp:
	${CC} -o ${prog_name}.exe ${prog_name}.cpp

run:
	./${prog_name}.exe


debug:
	gdb ./${prog_name}.exe
