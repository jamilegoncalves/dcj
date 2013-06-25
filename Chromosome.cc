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

/**
* Determinar número de adjacencias no cromossomo
*/

int Chromosome::numAdjacencies(Chromosome *chr)
{
    int numAdj;
    int numTel;

    int n = length();

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

