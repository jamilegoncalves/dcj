#include "AdjacencyGraph.h"
#include "Genome.h"
#include "Chromosome.h"

int main(int argc, char *argv[])
{
    		    Chromosome *a1 = new Chromosome("chrA1", true);

		    a1->genes.push_back(1);
		    a1->genes.push_back(19);
		    a1->genes.push_back(5);
		    a1->genes.push_back(-20);
		    a1->genes.push_back(3);

		    Chromosome *a2 = new Chromosome("chrA2", false);

		    a2->genes.push_back(4);
		    a2->genes.push_back(21);
		    a2->genes.push_back(2);

                    Chromosome *a3 = new Chromosome("chrA3", false);

                    a3->genes.push_back(26);
                    a3->genes.push_back(-27);

		    Genome *a = new Genome("Genome A");
		    a->chromosomes.push_back(a1);
		    a->chromosomes.push_back(a2);
                    a->chromosomes.push_back(a3);

                    Chromosome *b1 = new Chromosome("chrB1", true);

		    b1->genes.push_back(1);
		    b1->genes.push_back(23);
		    b1->genes.push_back(24);
		    b1->genes.push_back(2);
		    b1->genes.push_back(3);
		    b1->genes.push_back(25);
		    b1->genes.push_back(4);
		    b1->genes.push_back(5);

		    Genome *b = new Genome("Genome B");
		    b->chromosomes.push_back(b1);

		    AdjacencyGraph *ag = new AdjacencyGraph(a,b);

		    //std::cout << "Distancia pelo SortingByDCJ: " << ag->sortByDCJ() << std::endl;

                    //std::cout << "Distancia pela fórmula: " << ag->DCJdistance() << std::endl;

                    //std::cout << "Distancia pelo SortingByDCJRestrict: "
                            //<< ag->sortByRestrictedDCJ() << std::endl;
}
