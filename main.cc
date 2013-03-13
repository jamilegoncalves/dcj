#include "AdjacencyGraph.h"
#include "Genome.h"
#include "Chromosome.h"

int main(int argc, char *argv[])
{
	Chromosome *a1 = new Chromosome("chrA1", true);
        
        /*
        a1->genes.push_back(1);
        a1->genes.push_back(2);
        a1->genes.push_back(3);
        */
        
        
	a1->genes.push_back(1);
	a1->genes.push_back(3);
	a1->genes.push_back(-4);
        
        Chromosome *a2 = new Chromosome("chrA2", false);
	a2->genes.push_back(-2);
	a2->genes.push_back(-5);
        
        Chromosome *a3 = new Chromosome("chrA3", true);
	a3->genes.push_back(6);
	a3->genes.push_back(7);
        
        
	Genome *a = new Genome("Genome A");
        a->chromosomes.push_back(a1);
        a->chromosomes.push_back(a2);
        a->chromosomes.push_back(a3);
        
	Chromosome *b1 = new Chromosome("chrB1", false);
        
        /*
        b1->genes.push_back(1);
        b1->genes.push_back(-2);
        b1->genes.push_back(3);
        */
        
        
        b1->genes.push_back(-1);
	b1->genes.push_back(-2);

        Chromosome *b2 = new Chromosome("chrB2", true);
        b2->genes.push_back(3);
	b2->genes.push_back(4);
        
        Chromosome *b3 = new Chromosome("chrB3", true);
        b3->genes.push_back(5);

        Chromosome *b4 = new Chromosome("chrB4", false);
        b4->genes.push_back(-6);
	b4->genes.push_back(-7);
        
        
	Genome *b = new Genome("Genome B");
	b->chromosomes.push_back(b1);
        b->chromosomes.push_back(b2);
        b->chromosomes.push_back(b3);
        b->chromosomes.push_back(b4);

	std::cerr << *a << std::endl;
	std::cerr << *b << std::endl;
        
	AdjacencyGraph *ag = new AdjacencyGraph(a,b);

	std::cerr << *ag << std::endl;
        std::cout << "Distancia: " << ag->sortByDCJ() << std::endl;
}
