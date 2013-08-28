#include "AdjacencyGraph.h"
#include <stack>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Bergeron, Mixtacki, Stoye. A Unifying View
    // of Genome earrangements. 2006. Algorithm 1
    n = a->numGenes();
    numAdj = totalAdjacencies(a);

    adjA = adjB = NULL; // DEBUG

    idxEndOfAdjA = constructTables(a, adjA, locA)+1;
    idxEndOfAdjB = constructTables(b, adjB, locB)+1;

    //print(); // até aqui ok

//    paths(); // até aqui ok

//    capping(); // até aqui ok

//    print();

}

AdjacencyGraph::~AdjacencyGraph()
{
    delete adjA, adjB, locA, locB;
}

void AdjacencyGraph::print()
{
    printAdjacencies(std::cerr);
}

void printNumber(int num, std::ostream &os)
{
    if ( num > 0 )
    {
        os << num << "t";
    }
    else if ( num < 0 )
    {
        os << -num << "h";
    }
}

void AdjacencyGraph::printAdjacencies(std::ostream &os)
{
    if ( adjA != NULL )
    {
        os << "AdjA:" << std::endl;
        for(int i = 1; i < idxEndOfAdjA ; ++i)
        {
            printNumber(adjA[i].first, os);
            if (adjA[i].second != 0)
            {
                os << ",";
                printNumber(adjA[i].second, os);
            }
            os << " \t";

        }
        os << std::endl << std::endl;
    }

    if ( locA != NULL )
    {
        os << "LocA:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 16; ++i)
            os << "\t" << locA[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 16; ++i)
            os << "\t" << locA[i].tail;
        os << std::endl;
    }

    if ( adjB != NULL )
    {
        os << "AdjB:" << std::endl;
        for(int i = 1; i < idxEndOfAdjB; ++i)
        {
            printNumber(adjB[i].first, os);
            if (adjB[i].second != 0)
            {
                os << ",";
                printNumber(adjB[i].second, os);
            }
            os << " \t";
        }
        os << std::endl << std::endl;
    }

    if ( locB != NULL )
    {
        os << std::endl;
        os << "LocB:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 16; ++i)
            os << "\t" << locB[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 16; ++i)
            os << "\t" << locB[i].tail;
        os << std::endl;
    }
}

/*
void AdjacencyGraph::printAdjacencies(std::ostream &os)
{
    if ( adjA != NULL )
    {
        os << "AdjA:" << std::endl;
        os << "first:";
        for(int i = 1; i <= 22; ++i)
            os << "\t" << adjA[i].first;
        os << std::endl;

        os << "second:";
        for(int i = 1; i <= 22; ++i)
            os << "\t" << adjA[i].second;
        os << std::endl;
    }

    if ( locA != NULL )
    {
        os << "LocA:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 19; ++i)
            os << "\t" << locA[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 19; ++i)
            os << "\t" << locA[i].tail;
        os << std::endl;
    }

    if ( adjB != NULL )
    {
        os << std::endl;
        os << "AdjB:" << std::endl;
        os << "first:";
        for(int i = 1; i <= 22; ++i)
            os << "\t" << adjB[i].first;
        os << std::endl;

        os << "second:";
        for(int i = 1; i <= 22; ++i)
            os << "\t" << adjB[i].second;
        os << std::endl;
    }

    if ( locB != NULL )
    {
        os << std::endl;
        os << "LocB:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 19; ++i)
            os << "\t" << locB[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 19; ++i)
            os << "\t" << locB[i].tail;
        os << std::endl;
    }
}
*/

/**
 * Constroi tabelas: Adjacency e Location.
 * @returns Número de adjacências
 */
int AdjacencyGraph::constructTables(Genome *g, Adjacency *&adj, Location *&loc)
{
    int n = g->numGenes();
    int adjacencyTableSize = 3*n + 2;
    adj = new Adjacency[adjacencyTableSize];
    loc = new Location[2*n+1];
    int offset = 0;

    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        if(chr->isLinear() == true)
        {
            adj[offset + 1].first = (*chr)[1];
            adj[offset + 1].second = 0;
            adj[offset + 1].visited = false;
            adj[offset + chr->length() + 1].first = - (*chr)[chr->length()];
            adj[offset + chr->length() + 1].second = 0;
            adj[offset + chr->length() + 1].visited = false;

            for(int i = 2; i <= chr->length(); ++i)
            {
                adj[offset + i].first = - (*chr)[i-1];
                adj[offset + i].second = (*chr)[i];
                adj[offset + i].visited = false;
            }
        }
        else // if it's a circular chromosome
        {
            if(chr->length() == 1)
            {
                adj[offset + 1].first = (*chr)[1];
                adj[offset + 1].second = (*chr)[1];
                adj[offset + 1].visited = false;
            }
            else
            {
                adj[offset + 1].first = (*chr)[1];
                adj[offset + 1].second = - (*chr)[chr->length()];
                adj[offset + 1].visited = false;

                for(int i = 2; i <= chr->length(); ++i)
                {
                    adj[offset + i].first = - (*chr)[i-1];
                    adj[offset + i].second = (*chr)[i];
                    adj[offset + i].visited = false;
                }
            }
        }

        // Constroi tabela de locação
        for(int i = 1; i <= chr->length()+1; ++i)
        {
            int label = adj[offset + i].first;

            if(label > 0)
                loc[label].tail = offset + i;
            else if(label < 0)
                loc[-label].head = offset + i;

            label = adj[offset + i].second;

            if(label == 0)
                continue;
            else if(label > 0)
                loc[label].tail = offset + i;
            else if(label < 0)
                loc[-label].head = offset + i;
        }
        offset = offset + chr->numAdjacencies(chr);
    } // end for

    int numAdj = totalAdjacencies(g);

    for (int i = numAdj+1; i < adjacencyTableSize; ++i)
        adj[i].first = END_OF_TABLE;

    return numAdj;
}

/**
 * Implementation of Bergeron, Mixtacki & Stoye's greedy sorting
 * by DCJ (Algorithm 2).
 */
int AdjacencyGraph::sortByDCJ()
{
    Adjacency u, v, tempU, tempV;
    std::stack<int> vacancies;

    int dist = 0;

    // iterate over all adjacencies of genome B
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        // if it is an adjacency
        if(adjB[i].isAdjacency())
        {
            int idxU, idxV;
            int p = adjB[i].first, q = adjB[i].second;

            // let u be the element of genome A that contains p
            if(p > 0)
            {
                idxU = locA[p].tail;
                u = adjA[idxU];
            }
            else
            {
                idxU = locA[-p].head;
                u = adjA[idxU];
            }

            // let v be the element of genome A that contains q
            if(q > 0)
            {
                idxV = locA[q].tail;
                v = adjA[idxV];
            }
            else
            {
                idxV = locA[-q].head;
                v = adjA[idxV];
            }

            // if u != v then
            if( !u.equals(v) )
            {
                std::cout << "Cut: " << u.first << "," << u.second
                                                << std::endl;
                std::cout << "Cut: " << v.first << "," << v.second
                                                << std::endl;
                // replace u in A by {p,q}
                tempU.first = p;
                tempU.second = q;

                // replace v in A by (u\{p}) U (v\{q})
                tempV.first = u.setMinus(p);
                tempV.second = v.setMinus(q);

                if (tempV.first == 0)
                {
                    if(tempV.second == 0)
                        vacancies.push(idxV);
                    else
                    {
                        tempV.first = tempV.second;
                        tempV.second = 0;
                    }
                }
                std::cout << "Join: " << tempU.first << ","
                                            << tempU.second << std::endl;
                std::cout << "Join: " << tempV.first << ","
                                            << tempV.second << std::endl;

                // Altero a Tabela AdjA:
                adjA[idxU] = tempU;
                adjA[idxV] = tempV;

                // Altero a Tabela LocA:
                if(adjA[idxU].first > 0)
                    locA[adjA[idxU].first].tail = idxU;
                else
                    locA[-adjA[idxU].first].head = idxU;

                if(adjA[idxU].second > 0)
                    locA[adjA[idxU].second].tail = idxU;
                else
                    locA[-adjA[idxU].second].head = idxU;

                if(adjA[idxV].first > 0)
                    locA[adjA[idxV].first].tail = idxV;
                else
                    locA[-adjA[idxV].first].head = idxV;

                if(adjA[idxV].second > 0)
                    locA[adjA[idxV].second].tail = idxV;
                else
                    locA[-adjA[idxV].second].head = idxV;

                //print();
                ++dist;
                std::cout << "Distancia: " << dist << std::endl;
            }
        } // end if adjacency
    }// end for

    // iterate over all telomeres of genome B
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        // if it is a telomere in B
        if (adjB[i].isTelomere())
        {

            int idxU, idxV;
            int p = adjB[i].first;

            // let u be the element of genome A that contains p
            if(p > 0)
            {
                idxU = locA[p].tail;
                u = adjA[idxU];
            }
            else
            {
                idxU = locA[- p].head;
                u = adjA[idxU];
            }

            // if u is an adjacency
            if(u.isAdjacency())
            {
                std::cout << "Cut: " << u.first << "," << u.second
                                                << std::endl;

                // replace u in A by {p} ...
                tempU.first = p;
                tempU.second = 0;

                std::cout << "Join: " << tempU.first << ","
                                            << tempU.second << std::endl;

                // Altero Tabela AdjA:
                adjA[idxU] = tempU;

                // Altero Tabela LocA:
                if(tempU.first > 0)
                    locA[tempU.first].tail = idxU;
                else
                    locA[-tempU.first].head = idxU;

                // ... and (u\{p})
                tempV.first = u.setMinus(p);
                tempV.second = 0;

                std::cout << "Join: " << tempV.first << ","
                                            << tempV.second << std::endl;

                if(vacancies.empty())
                {
                    idxV = idxEndOfAdjA;
                    ++idxEndOfAdjA;
                }
                else
                {
                    idxV = vacancies.top();
                    vacancies.pop();
                }
                if(adjA[idxU].second > 0)
                    locA[adjA[idxU].second].tail = idxU;
                else
                    locA[-adjA[idxU].second].head = idxU;

                ++dist;

                //print();

                std::cout << "Distancia: " << dist << std::endl;
            } // end if u is an adjacency
        } // end if telomere
    }// end for

    return dist;
}

void AdjacencyGraph::prettyPrintA(std::ostream &os)
{
    bool visited[n+1];
    memset(&visited[0], 0, sizeof(visited));
    visited[0] = true;
    for (int i=1; i < n; ++i)
    {
        if (visited[i])
            continue;
        int k = i, previousk = 0;
        while (k != 0)
        {
            previousk = k;
            k = getPreviousinA(k);
        }
        for (int k = previousk; !visited[abs(k)]; k = getNextinA(k))
        {
            visited[abs(k)] = true;
            os << k << " ";
        }
        os << std::endl;
    }
    os << std::endl;
}

int AdjacencyGraph::sortByRestrictedDCJ()
{
    int distance=0;
    paths();
    capping();

    // Calcula e aplica a tabela de tradução:

    // Kovac: " Without loss of generality, we may assume that the markers
    // in chromosomes of Π are consecutive numbers ...
    // (otherwise renumber the markers). "

    std::map<int, int> translationTable, reverseTranslationTable;

    buildTranslationTable(adjB, locB, translationTable, reverseTranslationTable);

    std::map<int, int>::iterator entry_ = translationTable.begin();

    /*
    while ( entry_ != translationTable.end() )
    {
        std::pair<int, int> entry = *entry_;
        std::cout << entry.first << " -> " << entry.second << std::endl;
        ++entry_;
    }
     */

    // Renumber AdjB
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        adjB[i].first = translationTable[adjB[i].first];
        adjB[i].second = translationTable[adjB[i].second];
    }

    // Renumber AdjA
    for(int i = 1; i < idxEndOfAdjA; i++)
    {
        adjA[i].first = translationTable[adjA[i].first];
        adjA[i].second = translationTable[adjA[i].second];
    }

    rebuildLocation(adjA, locA, idxEndOfAdjA);
    rebuildLocation(adjB, locB, idxEndOfAdjB);

    print();

    // Procure na tabela de adjacência de A um telômero não visitado +x, 0:

    // Kovac: "We will be transforming Π into Γ gradually
    // ‘‘from left to right’’: once we have transformed the beginning of
    // a chromosome in Π to (ki , ki + 1, . . . , j ),
    // we extend it by moving j + 1 next to j. "

#ifdef DEBUG
    int oldDistance = -1;
#endif

    for (int j=1; j < n; j++)
    {
#ifdef DEBUG
        if (oldDistance != distance)
        {
            prettyPrintA(std::cerr);
            oldDistance = distance;
        }
#endif
        int idxj = locA[j].head;
        int idxjplus1 = locA[j+1].tail;

        // 1. if j+1 is next to j, we are done
        if ( !adjA[idxj].equals(adjB[locB[j].head]) )
        {
            // 2. if j, j+1 are in different chromosomes in A
            if(!sameChromosome(j,j+1))
            {
                // Translocation that brings j+1 close to j
                translocation(j);
                ++distance;
            }
            else
            {
                // 3. if j, j+1 have different orientations in A
                if(differentOrientationinA(j, j+1))
                {
                    // Reversal that brings j and j+1 together
                    reversal(j);
                    ++distance;
                }
                else
                {
                    // Otherwise, following Christie (1996), find the maker
                    // m with the highest number between j and j+1
                    int m = largestInABetween(j, j+1);

                    int idxm = locA[m].head;
                    int idxmplus1 = locA[m+1].tail;

                    // 4. If m+1 is on a different chromosome
                    if(!sameChromosome(m, m+1))
                    {
                        // Translocation to move m+1 next to m
                        translocation(m);
                        ++distance;

                        idxjplus1 = locA[j+1].tail; /// Alterei!!!!
                        // Translocation to move j+1 next to j
                        translocation(j);
                        print();
                        ++distance;
                    }
                    else
                    {
                        // Otherwise, the situation is j, ..., m, ..., j+1, ...
                        /// m+1

                        // 5. if m and m+1 have different orientation
                        if(differentOrientationinA(m, m+1))
                        {
                            // Reversal that brings m and m+1 together
                            reversal(m);
                            ++distance;
                            reversal(j);
                            ++distance;
                        }
                        else
                        {
                            // 6. Finally, m and m+1 have the same orientation
                            // Block interchange
                            if(m>0)
                            {
                                blockInterchange(j, m);
                                ++distance;
                                ++distance;
                            }
                            else
                            {
                                idxm = locA[m].tail;
                                idxmplus1 = locA[m+1].head;

                                blockInterchange(j, m);
                                ++distance;
                                ++distance;
                            }
                        }
                    }
                }
            }
        }
    }
    return distance;
}

bool AdjacencyGraph::differentOrientationinA(int left, int right)
{
    int i;

    for(i = left; abs(i) != abs(right); i = getNextinA(i))
    {
        assert(i != 0);
    }

    return ((left>0 && i<0)||(left<0 && i>0));
}

bool AdjacencyGraph::sameChromosome(int markerj, int markerk)
{
    int i;
    int j = 0;
    for(i = getNextinA(markerj); (i!= 0)&&(abs(i) != abs(markerk));
            i = getNextinA(i))
    {
        ++j;
    }

    if(i == 0)
        return false;
    else
        return true;

}

int AdjacencyGraph::getNextinA(int marker)
{
    int idx;

    if(marker > 0)
        idx = locA[marker].head;
    else
        idx = locA[-marker].tail;

    return adjA[idx].setMinus(-marker);
}

int AdjacencyGraph::getPreviousinA(int marker)
{
    int idx;

    if(marker > 0)
        idx = locA[marker].tail;
    else
        idx = locA[-marker].head;

    return -adjA[idx].setMinus(marker);
}

/**
 * Retorna o maior valor entre duas posições na tabela de adjacencias
 */
int AdjacencyGraph::largestInABetween(int left, int right)
{
    int maior = 0;
    for(int i = getNextinA(left); abs(i) != abs(right); i = getNextinA(i))
    {
        assert(i != 0);
        if( abs(i) >  abs(maior) )
            maior = i;
    }
    return maior;
}

/**
 * Reversal
 */

void AdjacencyGraph::reversal(int m)
{
    int idxj, idxk;
    int temp;

    idxj = locA[abs(m)].head;
    Adjacency *j = &adjA[idxj];

    idxk = locA[abs(m)+1].tail;
    Adjacency *k = &adjA[idxk];

    if( j->equals(adjA[idxj]) && k->equals(adjA[idxk]) )
    {
        if(adjA[idxj].first == m)
        {
            if(adjA[idxk].second == abs(m)+1)
            {
                temp = adjA[idxj].second;
                adjA[idxj].second = adjA[idxk].second;

                // Altero a Tabela LocA:
                if(adjA[idxk].second > 0)
                    locA[adjA[idxk].second].tail = idxj;
                else
                    locA[-adjA[idxk].second].head = idxj;

                adjA[idxk].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxk;
                else
                    locA[-temp].head = idxk;
            }
            else
            {
                temp = adjA[idxj].second;
                adjA[idxj].second = adjA[idxk].first;

                // Altero a Tabela LocA:
                if(adjA[idxk].first > 0)
                    locA[adjA[idxk].first].tail = idxj;
                else
                    locA[-adjA[idxk].first].head = idxj;

                adjA[idxk].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxk;
                else
                    locA[-temp].head = idxk;
            }
        }
        else
        {
            if(adjA[idxk].second == abs(m)+1)
            {
                temp = adjA[idxj].first;
                adjA[idxj].first = adjA[idxk].second;

                // Altero a Tabela LocA:
                if(adjA[idxk].second > 0)
                    locA[adjA[idxk].second].tail = idxj;
                else
                    locA[-adjA[idxk].second].head = idxj;

                adjA[idxk].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxk;
                else
                    locA[-temp].head = idxk;
            }
            else
            {
                temp = adjA[idxj].first;
                adjA[idxj].first = adjA[idxk].first;

                // Altero a Tabela LocA:
                if(adjA[idxk].first > 0)
                    locA[adjA[idxk].first].tail = idxj;
                else
                    locA[-adjA[idxk].first].head = idxj;

                adjA[idxk].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxk;
                else
                    locA[-temp].head = idxk;
            }
        }
    }
}

/**
 * Translocation
 */

void AdjacencyGraph::translocation(int m)
{
    int idxj, idxk;
    int temp;

    idxj = locA[abs(m)].head;
    Adjacency *j = &adjA[idxj];

    idxk = locA[abs(m)+1].tail;
    Adjacency *k = &adjA[idxk];

    if(adjA[idxj].first == m)
    {
        if(adjA[idxk].first == abs(m)+1)
        {
            temp = adjA[idxj].second;
            adjA[idxj].second = adjA[idxk].first;

            // Altero a Tabela LocA:
            if(adjA[idxk].first > 0)
                locA[adjA[idxk].first].tail = idxj;
            else
                locA[-adjA[idxk].first].head = idxj;

            adjA[idxk].first = temp;

            // Altero a Tabela LocA:
            if(temp > 0)
                locA[temp].tail = idxk;
            else
                locA[-temp].head = idxk;
        }
        else
        {
            temp = adjA[idxj].second;
            adjA[idxj].second = adjA[idxk].second;

            // Altero a Tabela LocA:
            if(adjA[idxk].second > 0)
                locA[adjA[idxk].second].tail = idxj;
            else
                locA[-adjA[idxk].second].head = idxj;

            adjA[idxk].second = temp;

            // Altero a Tabela LocA:
            if(temp > 0)
                locA[temp].tail = idxk;
            else
                locA[-temp].head = idxk;
        }
    }
    else
    {
        if(adjA[idxk].first == abs(m)+1)
        {
            temp = adjA[idxj].first;
            adjA[idxj].first = adjA[idxk].first;

            // Altero a Tabela LocA:
            if(adjA[idxk].first > 0)
                locA[adjA[idxk].first].tail = idxj;
            else
                locA[-adjA[idxk].first].head = idxj;

            adjA[idxk].first = temp;

            // Altero a Tabela LocA:
            if(temp > 0)
                locA[temp].tail = idxk;
            else
                locA[-temp].head = idxk;
        }
        else
        {
            temp = adjA[idxj].first;
            adjA[idxj].first = adjA[idxk].second;

            // Altero a Tabela LocA:
            if(adjA[idxk].second > 0)
                locA[adjA[idxk].second].tail = idxj;
            else
                locA[-adjA[idxk].second].head = idxj;

            adjA[idxk].second = temp;

            // Altero a Tabela LocA:
            if(temp > 0)
                locA[temp].tail = idxk;
            else
                locA[-temp].head = idxk;
        }
    }
}

/**
 * Block Interchange (criação do circular intermediário)
 */
void AdjacencyGraph::blockInterchange(int j, int m)
{
    int idxa, idxb, idxc, idxd;
    int temp;

    idxa = locA[abs(j)].head;
    Adjacency *a = &adjA[idxa];

    idxb = locA[abs(j)+1].tail;
    Adjacency *b = &adjA[idxb];

    idxc = locA[abs(m)].head;
    Adjacency *c = &adjA[idxc];

    idxd = locA[abs(m)+1].tail;
    Adjacency *d = &adjA[idxd];

    if(a->first > 0)
        idxa = locA[a->first].tail;
    else
        idxa = locA[-a->first].head;

    if(b->second > 0)
        idxb = locA[b->second].tail;
    else
        idxb = locA[-b->second].head;

    if(c->second > 0)
        idxc = locA[c->second].tail;
    else
        idxc = locA[-c->second].head;

    if(d->second > 0)
        idxd = locA[d->second].tail;
    else
        idxd = locA[-d->second].head;

    // trazer j+1 para j
    if( a->equals(adjA[idxa]) && c->equals(adjA[idxc]) )
    {
        if(adjA[idxa].first == j)
        {
            if(adjA[idxc].first == abs(j)+1)
            {
                temp = adjA[idxa].second;
                adjA[idxa].second = adjA[idxc].first;

                // Altero a Tabela LocA:
                if(adjA[idxc].first > 0)
                    locA[adjA[idxc].first].tail = idxa;
                else
                    locA[-adjA[idxc].first].head = idxa;

                adjA[idxc].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxc;
                else
                    locA[-temp].head = idxc;
            }
            else
            {
                temp = adjA[idxa].second;
                adjA[idxa].second = adjA[idxc].second;

                // Altero a Tabela LocA:
                if(adjA[idxc].second > 0)
                    locA[adjA[idxc].second].tail = idxa;
                else
                    locA[-adjA[idxc].second].head = idxa;

                adjA[idxc].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxc;
                else
                    locA[-temp].head = idxc;
            }
        }
        else
        {
            if(adjA[idxc].first == abs(j)+1)
            {
                temp = adjA[idxa].first;
                adjA[idxa].first = adjA[idxc].first;

                // Altero a Tabela LocA:
                if(adjA[idxc].first > 0)
                    locA[adjA[idxc].first].tail = idxc;
                else
                    locA[-adjA[idxc].first].head = idxc;

                adjA[idxc].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxa;
                else
                    locA[-temp].head = idxa;
            }
            else
            {
                temp = adjA[idxa].first;
                adjA[idxa].first = adjA[idxc].second;

                // Altero a Tabela LocA:
                if(adjA[idxc].second > 0)
                    locA[adjA[idxc].second].tail = idxc;
                else
                    locA[-adjA[idxc].second].head = idxc;

                adjA[idxc].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxa;
                else
                    locA[-temp].head = idxa;
            }
        }
    }

    // trazer m+1 para m
    if( b->equals(adjA[idxb]) && d->equals(adjA[idxd]) )
    {
        if(adjA[idxb].first == m)
        {
            if(adjA[idxd].first == abs(m)+1)
            {
                temp = adjA[idxb].second;
                adjA[idxb].second = adjA[idxd].first;

                // Altero a Tabela LocA:
                if(adjA[idxd].first > 0)
                    locA[adjA[idxd].first].tail = idxd;
                else
                    locA[-adjA[idxd].first].head = idxd;

                adjA[idxd].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxb;
                else
                    locA[-temp].head = idxb;
            }
            else
            {
                temp = adjA[idxb].second;
                adjA[idxb].second = adjA[idxd].second;

                // Altero a Tabela LocA:
                if(adjA[idxd].second > 0)
                    locA[adjA[idxd].second].tail = idxb;
                else
                    locA[-adjA[idxd].second].head = idxb;

                adjA[idxd].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxd;
                else
                    locA[-temp].head = idxd;
            }
        }
        else
        {
            if(adjA[idxd].first == abs(m)+1)
            {
                temp = adjA[idxb].first;
                adjA[idxb].first = adjA[idxd].first;

                // Altero a Tabela LocA:
                if(adjA[idxd].first > 0)
                    locA[adjA[idxd].first].tail = idxb;
                else
                    locA[-adjA[idxd].first].head = idxb;

                adjA[idxd].first = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxd;
                else
                    locA[-temp].head = idxd;
            }
            else
            {
                temp = adjA[idxb].first;
                adjA[idxb].first = adjA[idxd].second;

                // Altero a Tabela LocA:
                if(adjA[idxd].second > 0)
                    locA[adjA[idxd].second].tail = idxb;
                else
                    locA[-adjA[idxd].second].head = idxb;

                adjA[idxd].second = temp;

                // Altero a Tabela LocA:
                if(temp > 0)
                    locA[temp].tail = idxd;
                else
                    locA[-temp].head = idxd;
            }
        }
    }
}

/**
 * Reconstroi tabela de locação
 */

void AdjacencyGraph::rebuildLocation(Adjacency *adj, Location *loc, int n)
{
    for(int i=1; i < n; i++)
    {
        if(adj[i].first > 0)
            loc[adj[i].first].tail = i;
        else
            loc[-adj[i].first].head = i;

        if(adj[i].second > 0)
            loc[adj[i].second].tail = i;
        else if(adj[i].second < 0)
            loc[-adj[i].second].head = i;
    }
}

/**
 * Calcula a distância de DCJ pela fórmula
 * d = N - (C + I/2)
 */
int AdjacencyGraph::DCJdistance()
{
    paths();
    int i = oddPaths.size() + adjacencies.size();
    int c = cycles.size();

    int dist = n - (c + (i/2));

    return dist;
}

/**
 * Armazena:
 * Primeiro elemento de cada caminho ímpar na fila: oddpaths
 * Primeiro elemento de cada caminho par na fila: evenpaths
 */
void AdjacencyGraph::paths()
{
    for(int i=1; adjA[i].first != END_OF_TABLE; ++i)
        adjA[i].visited = false;

    for(int i=1; adjB[i].first != END_OF_TABLE; ++i)
        adjB[i].visited = false;

    // Para cada telomero do Genoma A
    for(int i=1; adjA[i].first != END_OF_TABLE; ++i)
    {
        if( adjA[i].isTelomere() && !adjA[i].visited )
        {
            int idxLast;
            int length = getLengthFromA(i, &idxLast);

            Path path = { adjA[i].first, idxLast };

            if( length%2 == 0 )
                evenPathsFromA.push(path);
            else {
                if(length == 1)
                    adjacencies.push(path);
                else
                    oddPaths.push(path);
            }
        }
    }

    // Para cada telomero do Genoma B
    for(int i=1; adjB[i].first != END_OF_TABLE; ++i)
    {
        if( adjB[i].isTelomere() && !adjB[i].visited )
        {
            int idxLast;
            int length = getLengthFromB(i, &idxLast);
            Path path = { adjB[i].first, idxLast};

            if( length%2 == 0 )
                evenPathsFromB.push(path);
            else
            {
                std::cerr << "Error: invalid adjacency table!" << std::endl;
                throw "Error: invalid adjacency table!";
            }

        }
    }

    // Se não é telomero
    for(int i=1; adjB[i].first != END_OF_TABLE; ++i)
    {
        if( !adjB[i].isTelomere() && !adjB[i].visited )
        {
            int length = getLengthFromB(i);
            cycles.push(adjB[i].first);
        }
    }

}

#define swapCond(a,b) if (a != 0 && b != 0) { temp = a; a = b; b = temp; }

int AdjacencyGraph::getLengthFromA(int i, int *idxLast)
{
    int length = 0;
    int j, k;
    int temp;

    if (idxLast != NULL) *idxLast = 0;

    do
    {
        adjA[i].visited = true;

        if(adjA[i].first > 0)
            j = locB[adjA[i].first].tail;
        else
            j = locB[-adjA[i].first].head;

        adjB[j].visited = true;

        length++;

        if(adjB[j].second == 0)
        {
            if (idxLast != NULL) *idxLast = j;
            break;
        }

        if(adjA[i].first == adjB[j].first)
        {
            if(adjB[j].second > 0)
                i = locA[adjB[j].second].tail;
            else
                i = locA[-adjB[j].second].head;

            if(adjB[j].second == adjA[i].first)
                swapCond(adjA[i].first, adjA[i].second);
        }
        else
        {
            if(adjB[j].first > 0)
                i = locA[adjB[j].first].tail;
            else
                i = locA[-adjB[j].first].head;

            if(adjB[j].first == adjA[i].first)
                swapCond(adjA[i].first, adjA[i].second);
        }

        length++;
    } while (adjA[i].second != 0 && !adjA[i].visited);

    if (idxLast != NULL && length % 2 == 0) *idxLast = i;

    return length;
}

int AdjacencyGraph::getLengthFromB(int i, int *idxLast)
{
    int length = 0;
    int j, k;
    int temp;

    if (idxLast != NULL) *idxLast = 0;

    do
    {
        adjB[i].visited = true;

        if(adjB[i].first > 0)
            j = locA[adjB[i].first].tail;
        else
            j = locA[-adjB[i].first].head;

        adjA[j].visited = true;

        length++;

        if(adjA[j].second == 0)
        {
            if (idxLast != NULL) *idxLast = j;
            break;
        }

        if(adjB[i].first == adjA[j].first)
        {
            if(adjA[j].second > 0)
                i = locB[adjA[j].second].tail;
            else
                i = locB[-adjA[j].second].head;

            if(adjA[j].second == adjB[i].first)
                swapCond(adjB[i].first, adjB[i].second);
        }
        else
        {
            if(adjA[j].first > 0)
                i = locB[adjA[j].first].tail;
            else
                i = locB[-adjA[j].first].head;

            if(adjA[j].first == adjB[i].first)
                swapCond(adjB[i].first, adjB[i].second);
        }

        length++;
    } while (adjB[i].second != 0 && !adjB[i].visited);

    if (idxLast != NULL && length % 2 == 0) {
        adjB[i].visited = true;
        *idxLast = i;
    }
    return length;
}

void AdjacencyGraph::buildTranslationTable(Adjacency *adj, Location *loc,
        std::map<int,int> &translationTable,
        std::map<int,int> &reverseTranslationTable)
{
    for(int j=1; adj[j].first != END_OF_TABLE; j++)
        adj[j].visited = false;

    int numGene = 0;

    for (int j=1; adj[j].first != END_OF_TABLE; j++)
    {
        if ( adj[j].visited == false && adj[j].second == 0 )
        {
            int last = 0;
            int k = j;
            int current;
            do	{
                current = adj[k].first;
                adj[k].visited = true;
                if ( current == -last )
                {
                    current = adj[k].second;
                }
                if ( current == 0 )
                    break;
                ++numGene;
                translationTable[current] = numGene;
                reverseTranslationTable[numGene] = current;
                translationTable[-current] = -numGene;
                reverseTranslationTable[-numGene] = -current;
                k = current > 0 ?
                        loc[current].head
                    :
                        loc[-current].tail; // vai para próximo gene no cromossomo
                last = current;
            }	while (true);
        }
    }
    translationTable[0] = 0;
    reverseTranslationTable[0] = 0;
}


/**
 * Insere os caps nos cromossomos lineares
 */
void AdjacencyGraph::capping()
{
    int c, p;
    int idxc1, idxc2, idxp1, idxp2;

    // Capping para caminhos ímpares:

    while(!oddPaths.empty())
    {
        c = oddPaths.front().start;

        // Adicionar:
        // (n+1)h no inicio do caminho
        // (n+1)t no final da tabela adjA

        if(c > 0)
            idxc1 = locA[c].tail;
        else
            idxc1 = locA[-c].head;

        adjA[idxc1].second = n+1;
        locA[n+1].tail = idxc1;
        adjA[idxEndOfAdjA].first = -(n+1);
        locA[n+1].head = idxEndOfAdjA;
        adjA[idxEndOfAdjA].second = 0;

        idxEndOfAdjA++;

        // Adicionar:
        // (n+1)h no final do caminho
        // (n+1)t no final da tabela adjB

        idxc2 = oddPaths.front().idxLast;

        adjB[idxc2].second = n+1;
        locB[n+1].tail = idxc2;
        adjB[idxEndOfAdjB].first = -(n+1);
        locB[n+1].head = idxEndOfAdjB;
        adjB[idxEndOfAdjB].second = 0;

        idxEndOfAdjB++;

        oddPaths.pop();

        cycles.push(c);
        n++;

    }

    // Capping para caminhos pares:

    while(!evenPathsFromA.empty() || !evenPathsFromB.empty())
    {
        if(!evenPathsFromA.empty())
        {
            p = evenPathsFromA.front().start;

            // Adicionar:
            // (n+1)h no inicio do caminho
            // (n+1)t no final da tabela adjA

            if(p > 0)
                idxp1 = locA[p].tail;
            else
                idxp1 = locA[-p].head;

            adjA[idxp1].second = n+1;
            locA[n+1].tail = idxp1;
            adjA[idxEndOfAdjA].first = -(n+1);
            locA[n+1].head = idxEndOfAdjA;
            adjA[idxEndOfAdjA].second = 0;

            idxEndOfAdjA++;

            // Adicionar:
            // (n+2)h no final do caminho
            // (n+2)t no final da tabela adjA

            idxp2 = evenPathsFromA.front().idxLast;

            adjA[idxp2].second = n+2;
            locA[n+2].tail = idxp2;
            adjA[idxEndOfAdjA].first = -(n+2);
            locA[n+2].head = idxEndOfAdjA;
            adjA[idxEndOfAdjA].second = 0;

            idxEndOfAdjA++;

            // Nova adjacência em adjB:
            adjB[idxEndOfAdjB].first = n+1;
            locB[n+1].tail = idxEndOfAdjB;
            adjB[idxEndOfAdjB].second = n+2;
            locB[n+2].tail = idxEndOfAdjB;

            idxEndOfAdjB++;

            // Adicionar:
            // (n+1)t no final da tabela adjB
            // (n+2)t no final da tabela adjB

            adjB[idxEndOfAdjB].first = -(n+1);
            locB[n+1].head = idxEndOfAdjB;
            adjB[idxEndOfAdjB].second = 0;

            idxEndOfAdjB++;

            adjB[idxEndOfAdjB].first = -(n+2);
            locB[n+2].head = idxEndOfAdjB;
            adjB[idxEndOfAdjB].first = 0;

            idxEndOfAdjB++;

            evenPathsFromA.pop();

            cycles.push(p);
            n = n+2;
        } // end if(!evenPathsFromA.empty())

        if(!evenPathsFromB.empty())
        {
            p = evenPathsFromB.front().start;

            // Adicionar:
            // (n+1)h no inicio do caminho
            // (n+1)t no final da tabela adjB

            if(p > 0)
                idxp1 = locB[p].tail;
            else
                idxp1 = locB[-p].head;

            adjB[idxp1].second = n+1;
            locB[n+1].tail = idxp1;
            adjB[idxEndOfAdjB].first = -(n+1);
            locB[n+1].head = idxEndOfAdjB;
            adjB[idxEndOfAdjB].second = 0;

            idxEndOfAdjB++;

            // Adicionar:
            // (n+2)h no final do caminho
            // (n+2)t no final da tabela adjB

            idxp2 = evenPathsFromB.front().idxLast;

            adjB[idxp2].second = n+2;
            locB[n+2].tail = idxp2;
            adjB[idxEndOfAdjB].first = -(n+2);
            locB[n+2].head = idxEndOfAdjB;
            adjB[idxEndOfAdjB].second = 0;

            idxEndOfAdjB++;

            // Nova adjacência:
            adjA[idxEndOfAdjA].first = n+1;
            locA[n+1].tail = idxEndOfAdjA;
            adjA[idxEndOfAdjA].second = n+2;
            locA[n+2].tail = idxEndOfAdjA;

            idxEndOfAdjA++;

            // Adicionar:
            // (n+1)t no final da tabela adjA
            // (n+2)t no final da tabela adjA


            adjA[idxEndOfAdjA].first = -(n+1);
            locA[n+1].head = idxEndOfAdjA;
            adjA[idxEndOfAdjA].second = 0;

            idxEndOfAdjA++;

            adjA[idxEndOfAdjA].first = -(n+2);
            locA[n+2].head = idxEndOfAdjA;
            adjA[idxEndOfAdjA].second = 0;

            idxEndOfAdjA++;

            evenPathsFromB.pop();

            cycles.push(p);
            n = n+2;
        } // end if(!evenPathsFromB.empty())
    }
}

/**
 * Determina número total de adjacencias no genoma.
 */
int AdjacencyGraph::totalAdjacencies(Genome *g)
{
    int totalAdj = 0;

    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        int n = chr->length();

        if(chr->isLinear() == true)
            totalAdj = totalAdj + n+1;
        else
            totalAdj = totalAdj + n;
    }
    return totalAdj;
}

bool Adjacency::equals(Adjacency &a)
{
    return((first == a.first)&&(second==a.second)
                            || (first==a.second)&&(second==a.first));
}

/* Returns true if it is an adjacency, false if it is a telomere.
 */
bool Adjacency::isAdjacency()
{
    if (first == 0)
        return false;
    if(second != 0)
        return true;
    else
        return false;
}

/* Returns true if it is a telomere, false if it is an adjacency.
 */
bool Adjacency::isTelomere()
{
    if (first == 0)
        return false;
    if(second == 0)
        return true;
    else
        return false;
}

/* returns this \ {x}  */
int Adjacency::setMinus(int x)
{
    if(first == x)
        return second;
    else
        return first;
}
