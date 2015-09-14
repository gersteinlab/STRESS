#CC = gcc -pg
CC = gcc
OBJS = GNCommunities.o List.o Flow.o Graph.o Tree.o
OBJS2 = SuboptimalPaths.o List.o Flow.o Graph.o Tree.o
OBJS3 = networkEvolution.o List.o Flow.o Graph.o alignment.o networkComparison.o Tree.o
OPTIONS = -c -O3 -Wall
OPTIONS = -c -Wall -ggdb -Wunused-but-set-variable
OOPTIONS = -Wall -O3 -o
OOPTIONS = -Wall -ggdb -o
LIBS = -lm
#OPTIONS = -c -O2 -m64  # for compiling on cardiff or cancun

all: gncommunities subopt

evolution: Evolution

gncommunities: $(OBJS)
	$(CC) $(OOPTIONS) gncommunities $(OBJS) $(LIBS)

subopt: $(OBJS2)
	$(CC) -o subopt $(OBJS2) $(LIBS)

treeTest: treeTest.o Tree.o
	$(CC) -o treeTest treeTest.o Tree.o $(LIBS)

Evolution: $(OBJS3)
	$(CC) -o Evolution $(OBJS3) $(LIBS)

List.o: List.c List.h
	$(CC) $(OPTIONS) List.c

Flow.o: Flow.c Flow.h
	$(CC) $(OPTIONS) Flow.c

Graph.o: Graph.c Graph.h
	$(CC) $(OPTIONS) Graph.c

GNCommunities.o: GNCommunities.c Graph.h
	$(CC) $(OPTIONS) GNCommunities.c

SuboptimalPaths.o: SuboptimalPaths.c Graph.h
	$(CC) $(OPTIONS) SuboptimalPaths.c

treeTest.o: treeTest.c
	$(CC) $(OPTIONS) treeTest.c

Tree.o: Tree.c Tree.h
	$(CC) $(OPTIONS) Tree.c

alignment.o: alignment.c alignment.h Graph.h
	$(CC) $(OPTIONS) alignment.c

networkComparison.o: networkComparison.c networkComparison.h alignment.h Graph.h List.h Community.h Flow.h
	$(CC) $(OPTIONS) networkComparison.c

networkEvolution.o: networkEvolution.c Graph.h alignment.h
	$(CC) $(OPTIONS) networkEvolution.c

clean:
	rm -f *.o *~ core
