#ifndef ADJACENCY_H
#define	ADJACENCY_H

#include "Genome.h"
#include <iostream>
#include <limits.h>
#include <utility>
#include <queue>
#include <map>
#include <vector>
#include <set>

class Adjacency
{
public:
    int first;
    int second;
    std::vector <int> label;
    std::pair<int,int> labelAdjWithFirst;
    std::pair<int,int> labelAdjWithSecond;

    bool visited;
    bool circularSingleton;

    Adjacency();
    bool isAdjacency();
    bool isTelomere();
    bool equals(Adjacency &a);

    /* returns this \ {x}  */
    int setMinus(int x);
};

#endif	/* ADJACENCY_H */