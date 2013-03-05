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