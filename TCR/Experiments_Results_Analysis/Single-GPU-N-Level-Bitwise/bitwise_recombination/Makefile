
objects = TNT.o TNT_gold.o

all: $(objects)
	nvcc -arch=sm_30 $(objects) -o app

TNT_gold.o: TNT_gold.cpp
	g++ -c TNT_gold.cpp
TNT.o: TNT.cu
	nvcc -c -I /uaopt/cuda/7.5.17/samples/common/inc -I. TNT.cu

clean:
	rm -f *.o app
