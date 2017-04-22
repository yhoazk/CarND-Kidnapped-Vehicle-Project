all: clean comp run

clean:
	rm  -rf ./*.exe

comp:
	g++ -o Motion_model.exe bayesianFilter.cpp main.cpp

run:
	./Motion_model.exe
