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

    int numCycles();
    int oddPaths();
    int sortByDCJ();

    int dcjDistance(Genome *a, Genome *b); // TODO: implementar usando formula

	/**
	 * Prints the representation of this AdjacencyGraph to a
	 * stream.
	 */
    friend std::ostream & operator<<( std::ostream &os,
	                                  const AdjacencyGraph& ag );

    private:

    Adjacency *adjA, *adjB;
    Location *locA, *locB;

    int sortByDCJ(int n, Adjacency *adjA, Adjacency *adjB,
            Location *locA, Location *locB);
                       // TODO: implementar Algorithm 2 e retornar dist

    void constructTables(Genome *g, Adjacency *&adj, Location *&loc);
    
    Genome *a;
    Genome *b;
    int n;
    Vertex *vertex;

};

#endif