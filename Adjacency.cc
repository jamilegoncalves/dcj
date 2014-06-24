#include "Adjacency.h"
#include "AdjacencyGraph.h"
#include <stack>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <set>
#include <bitset>
#include <algorithm>

Adjacency::Adjacency() :
first(0), second(0), visited(false), circularSingleton(false)
{
}

bool Adjacency::isAdjacency()
{
    if (first == 0)
        return false;
    if(second != 0)
        return true;
    else
        return false;
}

bool Adjacency::isTelomere()
{
    if (first == 0)
        return false;
    if(second == 0)
        return true;
    else
        return false;
}


bool Adjacency::equals(Adjacency &a)
{
    return((first == a.first)&&(second==a.second)
                            || (first==a.second)&&(second==a.first));
}

int Adjacency::setMinus(int x)
{
    if(first == x)
        return second;
    else
        return first;
}