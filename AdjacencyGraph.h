#ifndef ADJACENCYGRAPH_H
#define ADJACENCYGRAPH_H
#include "Genome.h"
#include <iostream>
#include <limits.h>

#define END_OF_TABLE INT_MAX

class Adjacency
{
public:
    int first;
    int second;
    bool equals(Adjacency& other);
    /* Returns true if it is an adjacency, false if it is a telomere.
     */
    bool isAdjacency();
    /* Returns true if it is a telomere, false if it is an adjacency.
     */
    bool isTelomere();

    /* returns this \ {x}  */
    int setMinus(int x);
};

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
    int idxEndOfAdjA;

    int sortByDCJ(int n, Adjacency *adjA, Adjacency *adjB,
            Location *locA, Location *locB);
                       // TODO: implementar Algorithm 2 e retornar dist

    /**
     * Constroi tabelas: Adjacency e Location.
     * @returns Número de adjacências
     */
    int constructTables(Genome *g, Adjacency *&adj, Location *&loc);

    Genome *a;
    Genome *b;
    int n, numAdj;
    Vertex *vertex;

};

#endif
