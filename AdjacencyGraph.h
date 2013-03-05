#ifndef ADJACENCYGRAPH_H
#define ADJACENCYGRAPH_H
#include "Genome.h"
#include <iostream>

typedef struct
{
    int first;
    int second;
} Adjacency;

typedef struct
{
    unsigned int head;
    unsigned int tail;
} Location;

typedef struct
{
    int vertex1;
    int vertex2;
} Vertex;

class AdjacencyGraph
{
    public:
    AdjacencyGraph(Genome *a, Genome *b);
    ~AdjacencyGraph();
    int totalAdjacencies(Genome *g);

    int sortByDCJ();
    
	/**
	 * Prints the representation of this AdjacencyGraph to a
	 * stream.
	 */
    friend std::ostream & operator<<( std::ostream &os,
	                                  const AdjacencyGraph& ag);

    void print();
    void printAdjacencies(std::ostream &os);

    private:

    Adjacency *adjA, *adjB;
    Location *locA, *locB;

    int sortByDCJ(int n, Adjacency *adjA, Adjacency *adjB,
            Location *locA, Location *locB);
                       // TODO: implementar Algorithm 2 e retornar dist

    void constructTables(Genome *g, Adjacency *&adj, Location *&loc);
    
    Genome *a;
    Genome *b;
    int n, numAdj;
    Vertex *vertex;

};

#endif
