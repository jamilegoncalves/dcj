#ifndef GENOME_H
#define GENOME_H
#include "Chromosome.h"
#include <vector>
#include <string>

class Genome
{
public:
	Genome(std::string name);
        int numGenes();

	std::vector<Chromosome *> chromosomes;
	std::string name;
};

#endif