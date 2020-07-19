#GXXFLAGS = -O3 -ffast-math -fprefetch-loop-arrays -mfpmath=sse -funswitch-loops -DNDEBUG -Wno-deprecated
GXXFLAGS = -O3 -Wno-deprecated -ggdb -m64 -DNDEBUG -march=k8 -ffast-math  -Wno-non-template-friend  -fprefetch-loop-arrays -mfpmath=sse -funroll-loops -funswitch-loops 

all:
	g++ -c -o spring.o source/spring.cpp $(GXXFLAGS)
	g++ -o spring $(GXXFLAGS) spring.o

linear:
	g++ -c -o spring_linear.o source/spring_linear.cpp $(GXXFLAGS)
	g++ -o spring_linear $(GXXFLAGS) spring_linear.o
clean:	
	rm -f *.o
	rm -f spring
