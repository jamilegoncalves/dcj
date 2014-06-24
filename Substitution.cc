#include "Substitution.h"
#include "AdjacencyGraph.h"
#include <iostream>
#include <queue>

Substitution::Substitution()
{

}

void Substitution::print(std::ostream &os) const
{
    os << "\n dcj substitution: (" << adj.first << "," << adj.second << ")";
    os << "\n";
    os << "\n";
}