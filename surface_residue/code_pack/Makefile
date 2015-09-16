CC = gcc -std=c99
OBJS = pdbFile.o modes.o transformations.o network.o surfaceProbe.o bindingSiteMeasures.o bindingSites.o

OPTIONS = -c -O3 -Wall -ggdb

LIBS = -lm  

all: bindingSites

bindingSites: $(OBJS)
	$(CC) -o bindingSites $(OBJS) $(LIBS)

bindingSites.o: bindingSites.c
	$(CC) $(OPTIONS) bindingSites.c $(LIBS)

pdbFile.o: pdbFile.c pdbFile.h
	$(CC) $(OPTIONS) pdbFile.c $(LIBS)

surfaceProbe.o: surfaceProbe.c surfaceProbe.h
	$(CC) $(OPTIONS) surfaceProbe.c $(LIBS)

network.o: network.c network.h
	$(CC) $(OPTIONS) network.c $(LIBS)

transformations.o: transformations.c transformations.h
	$(CC) $(OPTIONS) transformations.c $(LIBS)

modes.o: modes.c modes.h
	$(CC) $(OPTIONS) modes.c $(LIBS)

bindingSiteMeasures.o: bindingSiteMeasures.c bindingSiteMeasures.h
	$(CC) $(OPTIONS) bindingSiteMeasures.c $(LIBS)

clean:
	rm -f *.o core
