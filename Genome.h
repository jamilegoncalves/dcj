#ifndef GENOME_H
#define GENOME_H
#include "Chromosome.h"
#include <vector>
#include <string>
#include <iostream>

class Genome
{
public:
	Genome(std::string name);
        int numGenes();

        /**
         * Prints the representation of this Genome to a
         * stream.
         */
        friend std::ostream & operator<<( std::ostream &os,
                                          const Genome& g);

	std::vector<Chromosome *> chromosomes;
	std::string name;
};

#endif
