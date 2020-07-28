#GXXFLAGS = -O3 -ffast-math -fprefetch-loop-arrays -mfpmath=sse -funswitch-loops -DNDEBUG -Wno-deprecated
GXXFLAGS = -O3 -Wno-deprecated -ggdb -m64 -DNDEBUG -march=k8 -ffast-math  -Wno-non-template-friend  -fprefetch-loop-arrays -mfpmath=sse -funroll-loops -funswitch-loops -fpermissive -Wunused-result

all:
	g++ -c -o spring.o source/spring_interaction.cpp $(GXXFLAGS)
	g++ -o spring $(GXXFLAGS) spring.o
