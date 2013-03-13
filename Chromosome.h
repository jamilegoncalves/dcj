#ifndef CHROMOSOME_H
#define CHROMOSOME_H
#include <vector>
#include <string>
#include <iostream>

class Chromosome
{
public:
	Chromosome(std::string name, bool isLinear);

	std::vector<int> genes;
	std::string name;
	bool isLinear();
        int numAdjacencies(Chromosome *chr);

        /**
         * Prints the representation of this Genome to a
         * stream.
         */
        friend std::ostream & operator<<( std::ostream &os,
                                          const Chromosome& c);

        /**
         * Returns the i-th gene
         */
        int operator[](unsigned int i);

	/**
         * Returns the number of genes of this chromosome. */
	int length();
private:
        bool linear;
};

#endif
