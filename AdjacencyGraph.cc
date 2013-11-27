/* AdjacencyGraph.cc - Adjacency graph data structure
 * Copyright 2013 Jamile Gonçalves
 *
 * This file is part of DoubleCutAndJoin
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "AdjacencyGraph.h"
#include <stack>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <set>

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Braga, Machado, Ribeiro e Stoye. Genomic distance
    // under gene substitutions. 2011.

    n = a->numGenes();

    adjA = adjB = NULL;

    findLabels(a, b);
    
    idxEndOfAdjA = constructTables(a, labelsInA, adjA, locA, locLabelA)+1;
    idxEndOfAdjB = constructTables(b, labelsInB, adjB, locB, locLabelB)+1;

}

AdjacencyGraph::~AdjacencyGraph()
{
    delete adjA, adjB, locA, locB, locLabelA, locLabelB;
    if (labelsInA != NULL) delete labelsInA;
    if (labelsInB != NULL) delete labelsInB;
}

/**
 * Constroi tabelas: Adjacency e Location.
 * @returns Número de adjacências
 */
int AdjacencyGraph::constructTables(Genome *g, std::set<int> *labels,
        Adjacency *&adj,
        Location *&loc, LocationLabel *&locLabel)
{
    int n = g->numGenes();
    int adjacencyTableSize = 3*n + 2;

    adj = new Adjacency[adjacencyTableSize];
    loc = new Location[adjacencyTableSize];
    locLabel = new LocationLabel[adjacencyTableSize];
    int offset = 0;
    int numAdj = totalAdjacencies(g, labels);

    memset(locLabel, 0, adjacencyTableSize*sizeof(LocationLabel));

    std::vector<Chromosome*>::iterator cIterator;

    int k = 0;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        int l = labels[k].size();

        int x = 1;

        if(chr->isLinear() == true)
        {
            // Caso o cromossomo já comece com label
            while(labels[k].begin() != labels[k].end() &&
                    abs(*labels[k].begin()) == abs((*chr)[x]))
            {
                adj[offset+1].label.push_back((*chr)[x]);
                labels[k].erase((*chr)[x]);
                ++x;
            }

            if(x != chr->length()+1)
            {
                adj[offset+1].first = (*chr)[x];
                adj[offset+1].second = 0;

                x = 1;

                for(int j = x+1; j <= chr->numAdjacencies(chr, l); ++j)
                {
                    adj[offset + j].first = -(*chr)[x];
                    ++x;

                    while(labels[k].begin() != labels[k].end() &&
                            abs(*labels[k].begin()) == abs((*chr)[x]))
                    {
                        adj[offset + j].label.push_back((*chr)[x]);
                        labels[k].erase((*chr)[x]);
                        ++x;
                    }

                    adj[offset + j].second = (*chr)[x];
                } // End for
            } // End if is not Singleton
            else // if is LinearSingleton
            {
                adj[offset + 1].first = 0;
                adj[offset + 1].second = 0;
            }
        } // End if is Linear
        else // if it's a circular chromosome
        {
            if(chr->length() == 1)
            {
                if(!labels[k].empty())
                {
                    adj[offset + 1].first = 0;
                    adj[offset + 1].second = 0;
                    adj[offset + 1].label.push_back((*chr)[x]);
                    labels[k].erase((*chr)[x]);
                    circularSingleton.push_back(offset+1); // Armazena a posição
                                                      // do singleton circular
                    ++x;
                }
                else
                {
                    adj[offset + 1].first = (*chr)[1];
                    adj[offset + 1].second = -(*chr)[1];
                }
            }
            else
            {
                bool endLabel = false;
                // Caso o cromossomo já comece com label
                while(labels[k].begin() != labels[k].end() &&
                        abs(*labels[k].begin()) == abs((*chr)[x]))
                {
                    adj[offset+1].label.push_back((*chr)[x]);
                    labels[k].erase((*chr)[x]);
                    ++x;
                }

                if(x != chr->length()+1)
                {
                    adj[offset + 1].first = (*chr)[1];
                    // Caso que termine com label
                    while(labels[k].begin() != labels[k].end() &&
                            abs(*labels[k].begin()) == abs((*chr)[chr->length()]))
                    {
                        adj[offset+1].label.push_back((*chr)[chr->length()]);
                        labels[k].erase((*chr)[chr->length()]);
                        endLabel = true;
                        ++x;
                    }
                    if(endLabel == true)
                        adj[offset + 1].second = - (*chr)[x];
                    else
                        adj[offset + 1].second = - (*chr)[chr->length()];

                    x = 1;

                    for(int j = x+1; j <= chr->numAdjacencies(chr, l); ++j)
                    {
                        adj[offset + j].first = - (*chr)[x];
                        ++x;

                        while(labels[k].begin() != labels[k].end() &&
                                abs(*labels[k].begin()) == abs((*chr)[x]))
                        {
                            adj[offset + j].label.push_back((*chr)[x]);
                            labels[k].erase((*chr)[x]);
                            ++x;
                        }
                        
                        adj[offset + j].second = (*chr)[x];
                    }
                }
                else // if is a circular singleton
                {
                    adj[offset + 1].first = 0;
                    adj[offset + 1].second = 0;
                    circularSingleton.push_back(offset+1);
                }
            }
        } // end else is circular

        for (int i = numAdj+1; i < adjacencyTableSize; ++i)
            adj[i].first = END_OF_TABLE;

        // Constroi tabela de localização
        for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
        {
            if(adj[i].first > 0)
                loc[adj[i].first].tail = i;
            else
                loc[-adj[i].first].head = i;

            if(adj[i].second > 0)
                loc[adj[i].second].tail = i;
            else if(adj[i].second < 0)
                loc[-adj[i].second].head = i;

            if(!adj[i].label.empty())
            {
                for(std::vector<int>::iterator it = adj[i].label.begin();
                        it != adj[i].label.end(); it++)
                {
                    locLabel[abs(*it)].positionLabel = i;
                }
            }
        }
        ++k;
        offset = offset + chr->numAdjacencies(chr,l);
    } // end for each cromosome

    // Print
    std::cout<< "First: ";
    for(int i = 1; i <= 7; ++i)
        std::cout<< adj[i].first << ",";
    std::cout<< "\n";


    std::cout<< "Second: ";
    for(int i = 1; i <= 7; ++i)
        std::cout<< adj[i].second << ",";
    std::cout<< "\n";

    std::cout<< "Labels: ";
    for(int i = 1; i <= 7; ++i)
    {
        for (std::vector<int>::iterator it = adj[i].label.begin(); it != adj[i].label.end(); it++)
        {
            std::cout << *it << ",";
            // valor na posição apontada por it
        }
    }
    std::cout<< "\n";

    std::cout<< "Head: ";
    for(int i = 1; i <= 5; ++i)
        std::cout<< loc[i].head << ",";

    std::cout<< "\n";

    std::cout<< "tail: ";
    for(int i = 1; i <= 5; ++i)
        std::cout<< loc[i].tail << ",";
    std::cout<< "\n";

    std::cout<< "LocLabels: ";
    
    for (int i=19; i <= 27; ++i)
    {
        if(locLabel[i].positionLabel != 0)
            std::cout << locLabel[i].positionLabel << ",";
    }
    
    std::cout<< "\n";
    std::cout<< "\n";

}


bool Adjacency::isAdjacency()
{
    if (first == 0)
        return false;
    if(second != 0)
        return true;
    else
        return false;
}

bool Adjacency::isTelomere()
{
    if (first == 0)
        return false;
    if(second == 0)
        return true;
    else
        return false;
}

void AdjacencyGraph::findLabels(Genome *a, Genome *b)
{
    labels.clear();

    std::set<int> markersA, markersB;

    std::vector<Chromosome*>::iterator cIteratorA;
    std::vector<Chromosome*>::iterator cIteratorB;

    if (labelsInA != NULL) delete labelsInA;
    if (labelsInB != NULL) delete labelsInB;

    labelsInA = new std::set<int>[a->chromosomes.size()];
    labelsInB = new std::set<int>[b->chromosomes.size()];


    for(cIteratorA = a->chromosomes.begin();
            cIteratorA != a->chromosomes.end(); ++cIteratorA)
    {
        Chromosome *chrA = *cIteratorA;
        for(int i = 1; i <= chrA->length(); ++i)
            markersA.insert(abs((*chrA)[i]));
    }

    for(cIteratorB = b->chromosomes.begin();
            cIteratorB != b->chromosomes.end(); ++cIteratorB)
    {
        Chromosome *chrB = *cIteratorB;
        for(int i = 1; i <= chrB->length(); ++i)
            markersB.insert(abs((*chrB)[i]));
    }

    int k = 0;

    for(cIteratorA = a->chromosomes.begin();
            cIteratorA != a->chromosomes.end(); ++cIteratorA)
    {
        Chromosome *chrA = *cIteratorA;
        for(int i = 1; i <= chrA->length(); ++i)
        {
            if( markersB.find(abs((*chrA)[i])) == markersB.end() )
                labelsInA[k].insert(abs((*chrA)[i]));
        }
        ++k;
    }

    k = 0;

    for(cIteratorB = b->chromosomes.begin();
            cIteratorB != b->chromosomes.end(); ++cIteratorB)
    {
        Chromosome *chrB = *cIteratorB;
        for(int i = 1; i <= chrB->length(); ++i)
        {
            if( markersA.find(abs((*chrB)[i])) == markersA.end() )
                labelsInB[k].insert(abs((*chrB)[i]));
        }
        ++k;
    }
}

int AdjacencyGraph::totalAdjacencies(Genome *g, std::set<int> *labels)
{
    int totalAdj = 0;
    int k = 0;

    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        int n = chr->length();
        int l = labels[k].size();

        if(chr->isLinear() == true)
            totalAdj = totalAdj + n+1 - l;
        else
            totalAdj = totalAdj + n - l;
        
        if(n == l)
            totalAdj = totalAdj + 1;

        k++;
    }
    return totalAdj;
}

/*
bool Adjacency::equals(Adjacency &a)
{
    return((first == a.first)&&(second==a.second)
                            || (first==a.second)&&(second==a.first));
}

int Adjacency::setMinus(int x)
{
    if(first == x)
        return second;
    else
        return first;
}

*/