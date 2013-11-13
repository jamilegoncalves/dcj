#ifndef ADJACENCYGRAPH_H
#define ADJACENCYGRAPH_H
#include "Genome.h"
#include <iostream>
#include <limits.h>
#include <queue>
#include <map>
#include <vector>
#include <set>

#define END_OF_TABLE INT_MAX

class Adjacency
{
public:
    int first;
    int second;
    std::vector <int> label;

    //bool setCircularSingleton();
    //bool isCircularSingleton();
    //bool isLinearSingleton();
    bool isAdjacency();
    bool isTelomere();

private:
    bool circularSingleton;

/*
    bool equals(Adjacency& other);
    
    int setMinus(int x);
*/

};

typedef struct
{
    unsigned int head;
    unsigned int tail;
} Location;

typedef struct
{
    unsigned int label;
} LocationLabel;

class AdjacencyGraph
{
    public:
    AdjacencyGraph(Genome *a, Genome *b);
    ~AdjacencyGraph();

    private:


    Adjacency *adjA, *adjB;
    Location *locA, *locB;
    LocationLabel *locLabelA, *locLabelB;
    int idxEndOfAdjA, idxEndOfAdjB;

    /**
     * Constroi tabelas: Adjacency e Location.
     * @returns Número de adjacências
     */
    int constructTables(Genome *g,  std::set<int> *labels, Adjacency *&adj,
            Location *&loc, LocationLabel *&locLabel);

    void findLabels(Genome *a, Genome *b);

    Genome *a;
    Genome *b;
    std::set<int> labels;
    std::set<int> *labelsInA, *labelsInB;
    // Armazena a posição do Singleton Circular
    std::vector <int> circularSingleton;

    int n;

};

#endif
