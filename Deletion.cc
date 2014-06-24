#include "Deletion.h"
#include "AdjacencyGraph.h"
#include <iostream>

Deletion::Deletion()
{

}

void Deletion::print(std::ostream &os) const
{
    os << "\n dcj deletion: (" << adj.first << "," << adj.second << ")";
    os << "\n";
    os << "\n";
}