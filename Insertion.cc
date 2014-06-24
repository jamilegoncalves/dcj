#include "Insertion.h"
#include "AdjacencyGraph.h"
#include <iostream>

Insertion::Insertion()
{

}

void Insertion::print(std::ostream &os) const
{
    os << "\n dcj insertion: (" << adj.first << "," << adj.second << ")";
    os << "\n";
    os << "\n";
}