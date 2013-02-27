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
private:
        bool linear;
};

#endif