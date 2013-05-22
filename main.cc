#include "AdjacencyGraph.h"
#include "Genome.h"
#include "Chromosome.h"

int main(int argc, char *argv[])
{
		    Chromosome *a1 = new Chromosome("chrA1", true);

		    a1->genes.push_back(1);
		    a1->genes.push_back(-4);
		    a1->genes.push_back(2);
		    a1->genes.push_back(3);
		    a1->genes.push_back(5);
		    a1->genes.push_back(7);
		    a1->genes.push_back(6);
		    a1->genes.push_back(8);

		    Chromosome *a2 = new Chromosome("chrA2", true);
		    a2->genes.push_back(-17);
		    a2->genes.push_back(-14);
		    a2->genes.push_back(-15);
		    a2->genes.push_back(-13);
		    a2->genes.push_back(-11);
		    a2->genes.push_back(-12);
		    a2->genes.push_back(-10);
		    a2->genes.push_back(16);
		    a2->genes.push_back(-9);

		    Genome *a = new Genome("Genome A");
		    a->chromosomes.push_back(a1);
		    a->chromosomes.push_back(a2);

                    Chromosome *b1 = new Chromosome("chrB1", true);

		    b1->genes.push_back(2);
		    b1->genes.push_back(1);
		    b1->genes.push_back(3);
		    b1->genes.push_back(5);
		    b1->genes.push_back(4);

		    Chromosome *b2 = new Chromosome("chrB2", true);

                    b2->genes.push_back(6);
		    b2->genes.push_back(7);
		    b2->genes.push_back(-11);
		    b2->genes.push_back(-9);
		    b2->genes.push_back(-10);
		    b2->genes.push_back(-8);
		    b2->genes.push_back(12);
		    b2->genes.push_back(16);
		    b2->genes.push_back(17);

                    Chromosome *b3 = new Chromosome("chrB3", false);

                    b3->genes.push_back(15);
		    b3->genes.push_back(14);
		    b3->genes.push_back(-13);

		    Genome *b = new Genome("Genome B");
		    b->chromosomes.push_back(b1);
		    b->chromosomes.push_back(b2);
                    b->chromosomes.push_back(b3);

		    std::cerr << *a << std::endl;
		    std::cerr << *b << std::endl;

		    AdjacencyGraph *ag = new AdjacencyGraph(a,b);

		    std::cerr << *ag << std::endl;
		    std::cout << "Distancia: " << ag->sortByDCJ() << std::endl;

}