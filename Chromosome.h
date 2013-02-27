#ifndef CHROMOSOME_H
#define CHROMOSOME_H
#include <vector>
#include <string>

class Chromosome
{
public:
	Chromosome(std::string name, bool isLinear);

	std::vector<int> genes;
	std::string name;
	bool isLinear();
        int numAdjacencies(Chromosome *chr);

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
