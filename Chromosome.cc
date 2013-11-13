#include "Chromosome.h"

Chromosome::Chromosome(std::string name, bool isLinear)
{
    this->name = name;
    this->linear = isLinear;
    genes.push_back(0); // this is to make sure that
                        // the first index is 1
}

int Chromosome::operator[](unsigned int i)
{
    return genes[i];
}

bool Chromosome::isLinear()
{
    return linear;
}

int Chromosome::length()
{
    return genes.size()-1;
}

std::ostream & operator<<( std::ostream &os, const Chromosome& c)
{
    for (int i=1; i < c.genes.size(); ++i)
    {
        os << c.genes[i] << " ";
    }
    if ( c.linear )
        os << "| ";
    else
        os << ") ";

    return os;
}

/**
* Determinar número de adjacencias no cromossomo
*/

int Chromosome::numAdjacencies(Chromosome *chr, int numLabels)
{
    int numAdj;
    int numTel;

    int n = length();

    if(chr->isLinear() == true)
    {
        numAdj = n-numLabels-1;
        numTel = 2;
    }
    else
    {
        numAdj = n-numLabels;
        numTel = 0;
    }
    int offset = numAdj + numTel;
    return offset;
}
