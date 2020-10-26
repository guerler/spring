GXXFLAGS = -O3 -Wno-deprecated -Wno-gnu-static-float-init -ggdb -m64 -DNDEBUG -march=k8 -ffast-math -mfpmath=sse -funroll-loops -fpermissive -Wunused-result

all:
	rm -f *.o
	g++ -c -o spring.o source/spring.cpp $(GXXFLAGS)
	g++ -o spring $(GXXFLAGS) spring.o