#ifndef ADJACENCYGRAPH_H
#define ADJACENCYGRAPH_H
#include "Genome.h"
#include <iostream>
#include <limits.h>
#include <queue>
#include <map>

#define END_OF_TABLE INT_MAX

class Adjacency
{
public:
    int first;
    int second;
    bool visited;



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
    int start;
    int idxLast;
} Path;

class AdjacencyGraph
{
    public:
    AdjacencyGraph(Genome *a, Genome *b);
    ~AdjacencyGraph();
    int totalAdjacencies(Genome *g);

    int sortByDCJ();

    int DCJdistance();

    int sortByRestrictedDCJ();

/**
* Prints the representation of this AdjacencyGraph to a
* stream.
*/

    void print();
    void printAdjacencies(std::ostream &os);

    private:

    std::queue<int> cycles;
    std::queue<Path> oddPaths,evenPathsFromA, evenPathsFromB, adjacencies;

    Adjacency *adjA, *adjB;
    Location *locA, *locB;
    int idxEndOfAdjA, idxEndOfAdjB;
                       // TODO: implementar Algorithm 2 e retornar dist

    /**
     * Constroi tabelas: Adjacency e Location.
     * @returns Número de adjacências
     */
    int constructTables(Genome *g, Adjacency *&adj, Location *&loc);

    /**
     * Armazena:
     * Primeiro elemento de cada caminho ímpar na fila: oddpaths
     * Primeiro elemento de cada caminho par na fila: evenpaths
     */
    void paths();

    int getLengthFromA(int i, int *idxLast = NULL);

    int getLengthFromB(int i, int *idxLast = NULL);

    /**
     * Insere os caps nos cromossomos lineares
     */
    void capping();

    void buildTranslationTable(Adjacency *adj, Location *loc,
        std::map<int,int> &translationTable,
        std::map<int,int> &reverseTranslationTable);

    /**
     * Reconstroi tabela de locação
    */
    void rebuildLocation(Adjacency *adj, Location *loc, int n);

    Genome *a;
    Genome *b;
    int n, numAdj;

};

#endif
