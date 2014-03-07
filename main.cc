#include "AdjacencyGraph.h"
#include "Genome.h"
#include "Chromosome.h"

using namespace std;

int main(int argc, char *argv[])
{
    Genome *a;
    Genome *b;
    
    try {
        //a = new Genome("Genome A", "-1 19 | 2 20 3 | -7 -4  6 | 5 23 24 9 -11 12 | -8 -10 28 |");
        //b = new Genome("Genome B", "-1 -2 -3 21 -4 -7 22 -6 | -5 -9 25 11 -12 26 8 10 27 |");
        a = new Genome("Genome A", "1 19 5 -20 3 | 4 21 2 ) 26 -27 )");
        b = new Genome("Genome B", "1 23 24 2 3 25 4 5 |");
    } catch (const char *err) {
        std::cerr << err << std::endl;
        return 1;
    }
/*
    		    Chromosome *a1 = new Chromosome("chrA1", true);

		    a1->genes.push_back(-1);
		    a1->genes.push_back(19);

		    Chromosome *a2 = new Chromosome("chrA2", true);

		    a2->genes.push_back(2);
		    a2->genes.push_back(20);
		    a2->genes.push_back(3);

                    Chromosome *a3 = new Chromosome("chrA3", true);

                    a3->genes.push_back(-7);
                    a3->genes.push_back(-4);
                    a3->genes.push_back(6);

                    Chromosome *a4 = new Chromosome("chrA4", true);

                    a4->genes.push_back(5);
                    a4->genes.push_back(23);
                    a4->genes.push_back(24);
                    a4->genes.push_back(9);
                    a4->genes.push_back(-11);
                    a4->genes.push_back(12);

                    Chromosome *a5 = new Chromosome("chrA5", true);

                    a5->genes.push_back(-8);
                    a5->genes.push_back(-10);
                    a5->genes.push_back(28);

		    Genome *a = new Genome("Genome A");
		    a->chromosomes.push_back(a1);
		    a->chromosomes.push_back(a2);
                    a->chromosomes.push_back(a3);
                    a->chromosomes.push_back(a4);
                    a->chromosomes.push_back(a5);

                    Chromosome *b1 = new Chromosome("chrB1", true);

		    b1->genes.push_back(-1);
		    b1->genes.push_back(-2);
		    b1->genes.push_back(-3);
		    b1->genes.push_back(21);
		    b1->genes.push_back(-4);
		    b1->genes.push_back(-7);
		    b1->genes.push_back(22);
		    b1->genes.push_back(-6);

                    Chromosome *b2 = new Chromosome("chrB2", true);

		    b2->genes.push_back(-5);
		    b2->genes.push_back(-9);
		    b2->genes.push_back(25);
		    b2->genes.push_back(11);
		    b2->genes.push_back(-12);
		    b2->genes.push_back(26);
		    b2->genes.push_back(8);
		    b2->genes.push_back(10);
                    b2->genes.push_back(27);

		    Genome *b = new Genome("Genome B");
		    b->chromosomes.push_back(b1);
                    b->chromosomes.push_back(b2);
*/
		    AdjacencyGraph *ag = new AdjacencyGraph(a,b);

		    //std::cout << "Distancia pelo SortingByDCJ: " << ag->sortByDCJ() << std::endl;

                    //std::cout << "Distancia pela fórmula: " << ag->DCJdistance() << std::endl;

                    //std::cout << "Distancia pelo SortingByDCJRestrict: "
                            //<< ag->sortByRestrictedDCJ() << std::endl;
}
