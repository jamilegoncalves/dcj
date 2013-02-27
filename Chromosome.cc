#include "Chromosome.h"

Chromosome::Chromosome(std::string name, bool isLinear)
{
    this->name = name;
    this->linear = isLinear;
}

/**
 * Determinar número de adjacencias no cromossomo
 */

int Chromosome::numAdjacencies(Chromosome *chr)
{
    int numAdj;
    int numTel;
    
    int n = chr->genes.size();
        
    if(chr->isLinear() == true)
    {
        numAdj = n-1;
        numTel = 2;
    }
    else
    {
        numAdj = n;
        numTel = 0;
    }
    int offset = numAdj + numTel;
    return offset;
}

bool Chromosome::isLinear()
{
    return linear;
}