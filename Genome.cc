#include "Genome.h"

Genome::Genome(std::string name)
{
    this->name = name;
}

/**
 * Determinar número de genes no genoma
 */

int Genome::numGenes()
{
    int numGenes = 0;

    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = chromosomes.begin();
            cIterator != chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        numGenes = numGenes + chr->length();
    }
    return numGenes;
}

std::ostream & operator<<(std::ostream &os, const Genome& g)
{
    os << ">" << g.name << std::endl;

    for ( int i = 0; i < g.chromosomes.size(); ++i )
        os << *(g.chromosomes[i]);
    os << std::endl;

    return os;
}

