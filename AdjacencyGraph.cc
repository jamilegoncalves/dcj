#include "AdjacencyGraph.h"
#include <stack>

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Bergeron, Mixtacki, Stoye. A Unifying View
    // of Genome earrangements. 2006. Algorithm 1
    n = a->numGenes();
    numAdj = totalAdjacencies(a);
    vertex = new Vertex[2*n + 1];

    adjA = adjB = NULL; // DEBUG

    idxEndOfAdjA = constructTables(a, adjA, locA)+1;
    constructTables(b, adjB, locB);

    for(int i = 1; adjA[i].first != 0; ++i)
    {
        if(adjA[i].first > 0)
            vertex[i].vertex1 = locB[adjA[i].first].tail;
        else
            vertex[i].vertex1 = locB[- adjA[i].first].head;

        if(adjA[i].second == 0)
            vertex[i].vertex2 = 0;
        else if(adjA[i].second > 0)
            vertex[i].vertex2 = locB[adjA[i].second].tail;
        else
            vertex[i].vertex2 = locB[- adjA[i].second].head;
    }
    print();
}

AdjacencyGraph::~AdjacencyGraph()
{
    delete vertex;
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
        for(int i = 1; adjA[i].first != END_OF_TABLE ; ++i)
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
        os << "AdjB:" << std::endl;
        for(int i = 1; adjB[i].first != END_OF_TABLE ; ++i)
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
        for(int i = 1; i <= 19; ++i)
            os << "\t" << locB[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 19; ++i)
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

int AdjacencyGraph::sortByDCJ()
{
    return sortByDCJ(n, adjA, adjB, locA, locB);
}

/**
 * Constroi tabelas: Adjacency e Location.
 * @returns Número de adjacências
 */
int AdjacencyGraph::constructTables(Genome *g, Adjacency *&adj, Location *&loc)
{
    int n = g->numGenes();
    adj = new Adjacency[2*n + 2];
    loc = new Location[n+1];
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
            adj[offset + chr->length() + 1].first = - (*chr)[chr->length()];
            adj[offset + chr->length() + 1].second = 0;

            for(int i = 2; i <= chr->length(); ++i)
            {
                adj[offset + i].first = - (*chr)[i-1];
                adj[offset + i].second = (*chr)[i];
            }
        }
        else // if it's a circular chromosome
        {
            if(chr->length() == 1)
            {
                adj[offset + 1].first = (*chr)[1];
                adj[offset + 1].second = (*chr)[1];
            }
            else
            {
                adj[offset + 1].first = (*chr)[1];
                adj[offset + 1].second = - (*chr)[chr->length()];

                for(int i = 2; i <= chr->length(); ++i)
                {
                    adj[offset + i].first = - (*chr)[i-1];
                    adj[offset + i].second = (*chr)[i];
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

    for (int i = numAdj+1; i < 2*n+2; ++i)
        adj[i].first = END_OF_TABLE;

    return numAdj;
}

/**
 * Implementation of Bergeron, Mixtacki & Stoye's greedy sorting
 * by DCJ (Algorithm 2).
 */
int AdjacencyGraph::sortByDCJ(int n, Adjacency *adjA, Adjacency *adjB,
                                    Location *locA, Location *locB)
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

std::ostream & operator<<(std::ostream &os, const AdjacencyGraph& ag)
{
    for (int i = 1; i <= ag.numAdj; ++i)
        os << ag.vertex[i].vertex1 << "/"
           << ag.vertex[i].vertex2 << std::endl;

    return os;
}