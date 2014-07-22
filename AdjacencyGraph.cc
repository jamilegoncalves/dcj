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
#include "Adjacency.h"
#include "DoubleCutAndJoin.h"
#include "Deletion.h"
#include "Insertion.h"
#include "Substitution.h"
#include <stack>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <set>
#include <bitset>
#include <algorithm>

std::ostream & operator<<(std::ostream &os, const std::vector<int> &l);

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Braga, Machado, Ribeiro e Stoye. Genomic distance
    // under gene substitutions. 2011.

    n = a->numGenes();

    adjA = adjB = NULL;
    labelsInA = labelsInB = NULL;
    WhichGenome whereThis = undef;

    findLabels(a, b);

    whereThis = genomeA;
    idxEndOfAdjA = constructTables(a, labelsInA, adjA, adjAsize, locA, locAsize,
                                   locLabelA, locLabelAsize, whereThis)+1;
    whereThis = genomeB;
    idxEndOfAdjB = constructTables(b, labelsInB, adjB, adjBsize, locB, locBsize,
                                   locLabelB, locLabelBsize, whereThis)+1;
}

AdjacencyGraph::~AdjacencyGraph()
{
    delete adjA, adjB, locA, locB, locLabelA, locLabelB;
    if (labelsInA != NULL) delete labelsInA;
    if (labelsInB != NULL) delete labelsInB;
}

int AdjacencyGraph::maxGene(Genome *g)
{
    int max = -1;

    std::vector<Chromosome*>::iterator it;
    int k = 0;

    for(it = g->chromosomes.begin(); it != g->chromosomes.end(); ++it)
    {
        Chromosome *chr = *it;

        std::vector<int>::iterator itGene;

        for(itGene = chr->genes.begin(); itGene != chr->genes.end(); ++itGene)
        {
            if(abs(*itGene) > max)
                max = abs(*itGene);
        }
    }
    
    return max;
}

/**
 * Constroi tabelas: Adjacency e Location.
 * @returns Número de adjacências
 */
int AdjacencyGraph::constructTables(Genome *g, std::set<int> *labels,
        Adjacency *&adj, int &adjacencyTableSize,
        Location *&loc, int &locSize,
        LocationLabel *&locLabel, int &locLabelSize,
        WhichGenome whereThis)
{
    int n = g->numGenes();
    adjacencyTableSize = 3*n + 2;

    int maxG = maxGene(g);

    adj = new Adjacency[adjacencyTableSize]();
    locSize = maxG+1;
    loc = new Location[locSize];
    locLabelSize = maxG+1;
    locLabel = new LocationLabel[locLabelSize];
    int offset = 0;
    int numAdj = totalAdjacencies(g, labels);

    memset(locLabel, 0, (maxG+1)*sizeof(LocationLabel));
    memset(loc, 0, (maxG+1)*sizeof(Location));

    std::vector<Chromosome*>::iterator cIterator;

    int k = 0;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        int geneIdx;
        int numMarkers;

        if(chr->isLinear() == true)
        {
            geneIdx = 1;
            numMarkers = 0;

            // while this gene is a label...
            while (geneIdx <= chr->length() &&
                   labels[k].find(abs((*chr)[geneIdx])) != labels[k].end())
            {
                adj[offset+1].label.push_back((*chr)[geneIdx]);
                ++geneIdx;
            }
            
            // ALTEREI:
            //----------------------------------------------------------------------
            if(!adj[offset+1].label.empty())
            {
                // Preencho o labelAdjWithFirst:
                adj[offset+1].labelAdjWithFirst.first =
                        adj[offset+1].label.front();
                adj[offset+1].labelAdjWithFirst.second = 0;

                // Preencho o labelAdjWithSecond:
                adj[offset+1].labelAdjWithSecond.first =
                        -adj[offset+1].label.back();
                adj[offset+1].labelAdjWithSecond.second = (*chr)[geneIdx];
            }
            //----------------------------------------------------------------------

            if(geneIdx <= chr->length()) // if it's not a singleton
            {
                adj[offset+1].first = (*chr)[geneIdx];
                adj[offset+1].second = 0;
                ++numMarkers;

                while (geneIdx <= chr->length())
                {
                    adj[offset+numMarkers+1].first = -(*chr)[geneIdx];
                    adj[offset+numMarkers+1].second = 0;
                    
                    ++geneIdx;
                    // while this gene is a label...
                    while (geneIdx <= chr->length() &&
                        labels[k].find(abs((*chr)[geneIdx])) != labels[k].end()) {

                      adj[offset+numMarkers+1].label.push_back((*chr)[geneIdx]);
                      ++geneIdx;

                    }

                    // ALTEREI:
                    //----------------------------------------------------------------------
                    if(!adj[offset+numMarkers+1].label.empty())
                    {
                        // Preencho o labelAdjWithFirst:
                        adj[offset+numMarkers+1].labelAdjWithFirst.first =
                                -adj[offset+numMarkers+1].first;
                        adj[offset+numMarkers+1].labelAdjWithFirst.second =
                                adj[offset+numMarkers+1].label.front();

                        // Preencho o labelAdjWithSecond:
                        adj[offset+numMarkers+1].labelAdjWithSecond.first =
                                -adj[offset+numMarkers+1].second;
                        adj[offset+numMarkers+1].labelAdjWithSecond.second =
                                adj[offset+numMarkers+1].label.back();
                    }
                    //----------------------------------------------------------------------

                    ++numMarkers;
                    
                    if (geneIdx <= chr->length()) {
                        adj[offset+numMarkers].second = (*chr)[geneIdx];
                    }
                }

                offset += numMarkers;
            }
            else // if is LinearSingleton
            {
                adj[offset + 1].first = 0;
                adj[offset + 1].second = 0;
                
                if(whereThis == genomeA)
                    linearSingletonInA.push_back(offset+1);
                if(whereThis == genomeB)
                    linearSingletonInB.push_back(offset+1);

                ++offset;
            }
        } // End if is Linear
        else // if it's a circular chromosome
        {
            numMarkers = 0;

            int i = -1, j = -1; // i = index of the first marker
                                // j = index of the last marker
            // while this gene is a label...
            geneIdx = 1;

            while (geneIdx <= chr->length() &&
                   labels[k].find(abs((*chr)[geneIdx])) != labels[k].end())
            {
                adj[offset+1].label.push_back((*chr)[geneIdx]);
                ++geneIdx;
            }

            // ALTEREI:
            //----------------------------------------------------------------------

            if(!adj[offset+1].label.empty())
            {
                // Preencho o labelAdjWithFirst:
                adj[offset+1].labelAdjWithFirst.first =
                        adj[offset+1].second;
                adj[offset+1].labelAdjWithFirst.second =
                        adj[offset+1].label.front();

                // Preencho o labelAdjWithSecond:
                adj[offset+1].labelAdjWithSecond.first =
                        adj[offset+1].first;
                adj[offset+1].labelAdjWithSecond.second =
                        adj[offset+1].label.back();
            }
            //----------------------------------------------------------------------

            if(geneIdx <= chr->length()) {

                i = geneIdx;
                // iterate backwards
                geneIdx = chr->length();

                while (labels[k].find(abs((*chr)[geneIdx])) != labels[k].end())
                {
                    adj[offset+1].label.push_back((*chr)[geneIdx]);
                    --geneIdx;
                }

                j = geneIdx;

                // ALTEREI:
                //----------------------------------------------------------------------

                if(!adj[offset+1].label.empty())
                {
                    // Preencho o labelAdjWithSecond:
                    adj[offset+1].labelAdjWithSecond.first =
                            adj[offset+1].second;
                    adj[offset+1].labelAdjWithSecond.second =
                            -adj[offset+1].label.front();

                    // Preencho o labelAdjWithFirst:
                    adj[offset+1].labelAdjWithFirst.first =
                            adj[offset+1].first;
                    adj[offset+1].labelAdjWithFirst.second =
                            adj[offset+1].label.back();
                }
                //----------------------------------------------------------------------

                
            }
            if (i == -1) // if this chromosome has no markers
            {
                adj[offset+1].first = 0;
                adj[offset+1].second = 0;
                adj[offset+1].circularSingleton = true;

                if(whereThis == genomeA)
                    circularSingletonInA.push_back(geneIdx);
                if(whereThis == genomeB)
                    circularSingletonInB.push_back(geneIdx);

                ++offset;
            }
            else  // if this chromosome has at least one marker
            {
                geneIdx = i;
                adj[offset+1].first = (*chr)[i];
                adj[offset+1].second = -(*chr)[j];
                ++numMarkers;

                while(geneIdx < j)
                {
                    adj[offset + numMarkers + 1].first = -(*chr)[geneIdx];
                    ++geneIdx;

                    while (geneIdx <= j &&
                           labels[k].find(abs((*chr)[geneIdx])) != labels[k].end())
                    {
                        adj[offset + numMarkers + 1].label.push_back((*chr)[geneIdx]);
                        ++geneIdx;
                    }

                    adj[offset + numMarkers + 1].second = (*chr)[geneIdx];

                    // ALTEREI:
                    //----------------------------------------------------------------------

                    if(!adj[offset+numMarkers+1].label.empty())
                    {
                        // Preencho o labelAdjWithSecond:
                        adj[offset+numMarkers+1].labelAdjWithSecond.first =
                                adj[offset+numMarkers+1].second;
                        adj[offset+numMarkers+1].labelAdjWithSecond.second =
                                -adj[offset+numMarkers+1].label.front();

                        // Preencho o labelAdjWithFirst:
                        adj[offset+numMarkers+1].labelAdjWithFirst.first =
                                adj[offset+numMarkers+1].first;
                        adj[offset+numMarkers+1].labelAdjWithFirst.second =
                                adj[offset+numMarkers+1].label.back();
                    }
                    //----------------------------------------------------------------------

                    ++numMarkers;
                }
                offset += numMarkers;
            }
            // ALTEREI:
            //------------------------------------------------------------------
            if(!adj[1].label.empty())
            {
                // Preencho o labelAdjWithSecond:
                adj[1].labelAdjWithSecond.first =
                        adj[1].second;
                adj[1].labelAdjWithSecond.second =
                        adj[1].label.front();

                // Preencho o labelAdjWithFirst:
                adj[1].labelAdjWithFirst.first =
                        adj[1].first;
                adj[1].labelAdjWithFirst.second =
                        -adj[1].label.back();
            }
            //------------------------------------------------------------------

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
    } // end for each cromosome

    // Print
    std::cout<< "First: ";
    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
        std::cout<< adj[i].first << ",";
    std::cout<< "\n";


    std::cout<< "Second: ";
    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
        std::cout<< adj[i].second << ",";
    std::cout<< "\n";

    std::cout<< "Labels: ";
    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< adj[i].label << ",";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithFirst: ";
    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adj[i].labelAdjWithFirst.first << ","
                << adj[i].labelAdjWithFirst.second << "}" << ";";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithSecond: ";
    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adj[i].labelAdjWithSecond.first << ","
                << adj[i].labelAdjWithSecond.second << "}" << ";";
    }
    std::cout<< "\n";


    /*
    std::cout<< "Head: ";
    for(int i = 1; i <= 90; ++i)
        std::cout<< loc[i].head << ",";

    std::cout<< "\n";

    std::cout<< "tail: ";
    for(int i = 1; i <= 90; ++i)
        std::cout<< loc[i].tail << ",";
    std::cout<< "\n";
/*
    std::cout<< "LocLabels: ";
    
    for (int i=-20; i <= 25; ++i)
    {
        if(locLabel[i].positionLabel != 0)
            std::cout << locLabel[i].positionLabel << ",";
    }
*/

    std::cout<< "\n";
    std::cout<< "\n";

    return offset;
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

            Path path = { genomeA, adjA[i].first, idxLast };

            //adjA[idxLast].visited = true;

            if( length%2 == 0 )
                evenPathsFromA.push_back(path);
            else {
                if(length == 1)
                    pathsLength1.push_back(path);
                else
                    oddPaths.push_back(path);
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

            Path path = { genomeB, adjB[i].first, idxLast};
            //adjB[idxLast].visited = true;

            if( length%2 == 0 )
                evenPathsFromB.push_back(path);
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
            std::pair<WhichGenome,int> p(genomeB, adjB[i].first);
            cycles.push_back(p);
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

    if (idxLast != NULL && length % 2 == 0)
    {
        adjA[i].visited = true;
        *idxLast = i;
    }

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

Genome *AdjacencyGraph::adjacencyTableToGenome(Adjacency *adj, Location *loc)
{
    Genome *g = new Genome();
    std::queue< Adjacency > chromosomeRepresentatives;

    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
        adj[i].visited = false;

    int idxChr;
    int idxAdj = 1;
    bool isLinear;

    for(int i = 1; adj[i].first != END_OF_TABLE; ++i)
    {
        if(!adj[i].visited)
        {
            idxChr = i;
            int currentGene = adj[idxChr].first;

            if (currentGene == 0) // if it's a singleton
            {
                adj[i].visited = true;
                if (!adj[i].label.empty())
                    chromosomeRepresentatives.push(adj[i]);
                
                continue;
            }

            do
            {
                adj[idxChr].visited = true;

                int nextGene = adj[idxChr].setMinus(currentGene);

                if (nextGene > 0)
                    idxChr = loc[nextGene].head;
                else if (nextGene < 0)
                    idxChr = loc[-nextGene].tail;

                currentGene = -nextGene;
            }
            while( (!adj[idxChr].visited) && (currentGene != 0) );

            chromosomeRepresentatives.push(adj[idxChr]);

            if (currentGene == 0) // se for chr linear
            {                     // procura o outro telomero
                idxChr = i;
                currentGene = adj[idxChr].second;
                do
                {
                    currentGene = adj[idxChr].setMinus(-currentGene);

                    if(currentGene > 0)
                        idxChr = loc[currentGene].head;
                    else if(currentGene < 0)
                        idxChr = loc[-currentGene].tail;

                    adj[idxChr].visited = true;
                }
                while (currentGene != 0);
            }
        }
    }
    
    std::cout<< "Número de representantes: ";
    std::cout<<chromosomeRepresentatives.size();
    std::cout<< "\n";

    while (!chromosomeRepresentatives.empty())
    {
        Adjacency rep = chromosomeRepresentatives.front();
        chromosomeRepresentatives.pop();

        if(rep.first == 0) // If it's a singleton
        {
            if(rep.circularSingleton == true) // If it's a linear singleton
                isLinear = true;
            if(rep.circularSingleton == false) // If it's a circular singleton
                isLinear = false;
        }
        else // If it's not a singleton
        {
            if(rep.second == 0)
                isLinear = true;
            else
                isLinear = false;
        }

        Chromosome *chr = new Chromosome("chrC", isLinear);

        if(isLinear) // if it's a linear chromosome
        {
            if( (rep.first == 0)&&(!rep.label.empty()) ) // If it's a singleton
            {
                for(std::vector<int>::iterator it = rep.label.begin();
                            it != rep.label.end(); it++)
                {
                    chr->genes.push_back(*it);
                }
            }
            else // If it's not a singleton
            {
                int currentGene = rep.first;

                chr->genes.push_back(currentGene);

                if(!rep.label.empty())
                {
                    for(std::vector<int>::iterator it = rep.label.begin();
                            it != rep.label.end(); it++)
                    {
                        chr->genes.push_back(*it);
                    }
                }

                // next adjacency
                if(currentGene > 0)
                    idxAdj = loc[currentGene].head;
                else
                    idxAdj = loc[-currentGene].tail;

                int nextGene = adj[idxAdj].setMinus(-currentGene);

                do{
                    currentGene = nextGene;

                    if(!adj[idxAdj].label.empty())
                    {
                        for(std::vector<int>::iterator it = adj[idxAdj].label.begin();
                                it != adj[idxAdj].label.end(); it++)
                        {
                            chr->genes.push_back(*it);
                        }
                    }

                    chr->genes.push_back(currentGene);

                    // next adjacency
                    if(currentGene > 0)
                        idxAdj = loc[currentGene].head;
                    else
                        idxAdj = loc[-currentGene].tail;

                    nextGene = adj[idxAdj].setMinus(-currentGene);

                }while(nextGene != 0);

            }
            g->chromosomes.push_back(chr);
        }
        else // if it's a circular chromosome
        {
            if( (rep.first == 0)&&(!rep.label.empty()) ) // If it's a singleton
            {
                for(std::vector<int>::iterator it = rep.label.begin();
                            it != rep.label.end(); it++)
                {
                    chr->genes.push_back(*it);
                }
            }
            else // If it's not a singleton
            {
                int currentGene = rep.first;

                if(!rep.label.empty())
                {
                    for(std::vector<int>::iterator it = rep.label.begin();
                            it != rep.label.end(); it++)
                    {
                        chr->genes.push_back(*it);
                    }
                }

                chr->genes.push_back(currentGene);

                // next adjacency
                if(currentGene > 0)
                    idxAdj = loc[currentGene].head;
                else
                    idxAdj = loc[-currentGene].tail;

                int nextGene = adj[idxAdj].setMinus(-currentGene);

                do{

                    currentGene = nextGene;

                    if(!adj[idxAdj].label.empty())
                    {
                        for(std::vector<int>::iterator it = adj[idxAdj].label.begin();
                                it != adj[idxAdj].label.end(); it++)
                        {
                            chr->genes.push_back(*it);
                        }
                    }
                    
                    chr->genes.push_back(currentGene);

                    // next adjacency
                    if(currentGene > 0)
                        idxAdj = loc[currentGene].head;
                    else
                        idxAdj = loc[-currentGene].tail;

                    nextGene = adj[idxAdj].setMinus(-currentGene);

                }while(nextGene != rep.first);
            }
            g->chromosomes.push_back(chr);
        }
    }
    return g;
}

void AdjacencyGraph::printGenome(Genome *g)
{
    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = g->chromosomes.begin();
            cIterator != g->chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;
        std::cout<< *chr << std::endl;
    }
}

void alterTable(int idxAdj, int idxAdjB, Adjacency* &adj, int &idxEndOfAdj, int &adjSize,
                Location* &loc, int &locSize, Adjacency* &adjB)
{
        adj[idxAdj].first =  adjB[idxAdjB].labelAdjWithFirst.first;
        adj[idxAdj].second =  adjB[idxAdjB].labelAdjWithFirst.second;

        // Altera tabela de localização
        if(adj[idxAdj].first > 0)
            loc[adj[idxAdj].first].tail = idxAdj;
        else
            loc[-adj[idxAdj].first].head = idxAdj;

        if(adj[idxAdj].second > 0)
            loc[adj[idxAdj].second].tail = idxAdj;
        else
            loc[-adj[idxAdj].second].head = idxAdj;

        if (idxEndOfAdj >= adjSize) {
            // acabou o espaço, realocar
            Adjacency *newAdj = new Adjacency[2*adjSize];
            memcpy(newAdj, adj, sizeof(Adjacency) * adjSize);
            adjSize = 2*adjSize;
            delete adj;
            adj = newAdj;
        }

        adj[idxEndOfAdj].first = adjB[idxAdjB].labelAdjWithSecond.first;
        adj[idxEndOfAdj].second = adjB[idxAdjB].labelAdjWithSecond.second;

        // Altera tabela de localização
        if(adj[idxEndOfAdj].second > 0)
            loc[adj[idxEndOfAdj].second].tail = idxEndOfAdj;
        else
            loc[-adj[idxEndOfAdj].second].head = idxEndOfAdj;

        if(adj[idxEndOfAdj].first > 0)
            loc[adj[idxEndOfAdj].first].tail = idxEndOfAdj;
        else
            loc[-adj[idxEndOfAdj].first].head = idxEndOfAdj;
        
        ++idxEndOfAdj;
}

/*
void alterTable(int idxAdj, int label, Adjacency* &adj, int &idxEndOfAdj, int &adjSize,
                Location* &loc, int &locSize, Adjacency* &adjB) {
    alterTable(idxAdj, label, adj, idxEndOfAdj, adjSize, loc, locSize);

    adj[idxAdj].first =  adjB[idxAdj].labelAdjWithFirst.first;
    adj[idxAdj].second =  adjB[idxAdj].labelAdjWithFirst.second;

    adj[idxEndOfAdj].first = adjB[idxAdj].labelAdjWithSecond.first;
    adj[idxEndOfAdj].second = adjB[idxAdj].labelAdjWithSecond.second;
}
*/

/**
 * Cria uma nova adjacencia nas Tabelas adjA e adjB para as operações:
 * inserção e substituição
 */
void AdjacencyGraph::newAdjacency(int idxAdjA, int idxAdjB)
{
    int temp;

    for(std::vector<int>::iterator it = adjB[idxAdjB].label.begin();
                    it != adjB[idxAdjB].label.end(); it++)
    {
        alterTable(idxAdjA, idxAdjB, adjA, idxEndOfAdjA, adjAsize, locA, locAsize, adjB);
    }
    adjA[idxAdjA].label.clear();

    // Imprime tabela de adjacencias:
    // Print
    std::cout<< "\n";
    std::cout<< "First: ";
    for(int i = 1; adjA[i].first != END_OF_TABLE; ++i)
        std::cout<< adjA[i].first << ",";
    std::cout<< "\n";


    std::cout<< "Second: ";
    for(int i = 1; adjA[i].first != END_OF_TABLE; ++i)
        std::cout<< adjA[i].second << ",";
    std::cout<< "\n";

    std::cout<< "Labels: ";
    for(int i = 1; adjA[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< adjA[i].label << ",";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithFirst: ";
    for(int i = 1; adjA[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adjA[i].labelAdjWithFirst.first << ","
                << adjA[i].labelAdjWithFirst.second << "}" << ";";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithSecond: ";
    for(int i = 1; adjA[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adjA[i].labelAdjWithSecond.first << ","
                << adjA[i].labelAdjWithSecond.second << "}" << ";";
    }
    std::cout<< "\n";


    for(std::vector<int>::iterator it = adjB[idxAdjB].label.begin();
                    it != adjB[idxAdjB].label.end(); it++)
    {
        alterTable(idxAdjB, idxAdjB, adjB, idxEndOfAdjB, adjBsize, locB, locBsize, adjB);
    }
    adjB[idxAdjB].label.clear();

    // Imprime Tabela de Adjacencias
    // Print

    std::cout<< "\n";
    std::cout<< "First: ";
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
        std::cout<< adjB[i].first << ",";
    std::cout<< "\n";

    std::cout<< "Second: ";
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
        std::cout<< adjB[i].second << ",";
    std::cout<< "\n";

    std::cout<< "Labels: ";
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< adjB[i].label << ",";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithFirst: ";
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adjB[i].labelAdjWithFirst.first << ","
                << adjB[i].labelAdjWithFirst.second << "}" << ";";
    }
    std::cout<< "\n";

    std::cout<< "labelAdjWithSecond: ";
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        std::cout<< "{" << adjB[i].labelAdjWithSecond.first << ","
                << adjB[i].labelAdjWithSecond.second << "}" << ";";
    }
    std::cout<< "\n";

}

int AdjacencyGraph::sortByDCJsubst(std::queue<Genome *> &steps,
        std::queue< Rearrangement > &dcjs)
{
    Adjacency u, v, tempU, tempV;
    std::stack<int> vacancies;

    int dist = 0;

    // Print
    Genome *g = adjacencyTableToGenome(adjA, locA);
    //printGenome(g);

    steps.push(g);
    
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
                if(!adjB[i].label.empty())
                {
                    // Caso 1
                    if( (v.second == q) && (u.first == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempU.label.push_back(*it);
                        }
                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempU.label.push_back(*it);                            
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = u.setMinus(p);
                        tempV.second = v.setMinus(q);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);

                    }

                    // Caso 2
                    if( (v.first == q) && (u.first == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempU.label.push_back(*it);
                        }
                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempU.label.push_back(-(*it));
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = u.setMinus(p);
                        tempV.second = v.setMinus(q);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }

                    // Caso 3
                    if( (v.second == q) && (u.second == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempU.label.push_back(-(*it));
                        }
                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempU.label.push_back(*it);
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = u.setMinus(p);
                        tempV.second = v.setMinus(q);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }

                    // Caso 4
                    if( (v.first == q) && (u.second == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempU.label.push_back(-(*it));
                        }
                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempU.label.push_back(-(*it));                            
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = u.setMinus(p);
                        tempV.second = v.setMinus(q);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }
                }
                else
                {
                    // Caso 5
                    if( (v.second == q) && (u.first == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempU.label.push_back(*it);
                        }
                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempV.label.push_back(*it);
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = v.setMinus(q);
                        tempV.second = u.setMinus(p);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }

                    // Caso 6
                    if( (v.first == q) && (u.first == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempV.label.push_back(-(*it));                            
                        }
                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempV.label.push_back(*it);
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = v.setMinus(q);
                        tempV.second = u.setMinus(p);
                        
                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }

                    // Caso 7
                    if( (v.second == q) && (u.second == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempV.label.push_back(*it);
                        }
                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempV.label.push_back(-(*it));                            
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = v.setMinus(q);
                        tempV.second = u.setMinus(p);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }

                    // Caso 8
                    if( (v.first == q) && (u.second == p) )
                    {
                        DoubleCutAndJoin dcj;

                        dcj.cut[0] = u;
                        dcj.cut[1] = v;

                        // replace u in A by {p,q}
                        tempU.first = p;
                        tempU.label.clear();
                        tempU.second = q;

                        for(std::vector<int>::iterator it = v.label.begin();
                                it != v.label.end(); it++)
                        {
                            tempV.label.push_back(-(*it));                            
                        }
                        for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                        {
                            tempV.label.push_back(-(*it));
                        }

                        // replace v in A by (u\{p}) U (v\{q})
                        tempV.first = v.setMinus(q);
                        tempV.second = u.setMinus(p);

                        dcj.join[0] = tempU;
                        dcj.join[1] = tempV;

                        dcj.print(std::cerr);

                        dcjs.push(dcj);
                    }
                }

                if (tempV.first == 0)
                {
                    if( (tempV.second == 0) && (!tempV.label.empty()) )
                        vacancies.push(idxV);
                    else
                    {
                        tempV.first = tempV.second;
                        tempV.second = 0;
                    }
                }

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

                if(!adjA[idxU].label.empty())
                {
                    for(std::vector<int>::iterator it = adjA[idxU].label.begin();
                            it != adjA[idxU].label.end(); it++)
                    {
                        locLabelA[abs(*it)].positionLabel = idxU;
                    }
                }

                if(!adjA[idxV].label.empty())
                {
                    for(std::vector<int>::iterator it = adjA[idxV].label.begin();
                            it != adjA[idxV].label.end(); it++)
                    {
                        locLabelA[abs(*it)].positionLabel = idxV;
                    }
                }

                // Print Genome
                Genome *g = adjacencyTableToGenome(adjA, locA);
                //printGenome(g);

                steps.push(g);

                ++dist;
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
                if(!adjB[i].label.empty())
                {
                    DoubleCutAndJoin dcj;

                    dcj.cut[0] = u;
                    dcj.cut[1] = v;

                    // replace u in A by {p} ...
                    tempU.first = p;
                    tempU.label.clear();
                    tempU.second = 0;

                    for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                    {
                        tempU.label.push_back(*it);
                    }

                    // ... and (u\{p})
                    tempV.first = u.setMinus(p);
                    tempV.second = 0;

                    dcj.join[0] = tempU;
                    dcj.join[1] = tempV;

                    dcjs.push(dcj);
                }
                else
                {
                    DoubleCutAndJoin dcj;

                    dcj.cut[0] = u;
                    dcj.cut[1] = v;

                    // replace u in A by {p} ...
                    tempU.first = p;
                    tempU.label.clear();
                    tempU.second = 0;

                    for(std::vector<int>::iterator it = u.label.begin();
                                it != u.label.end(); it++)
                    {
                        tempV.label.push_back(*it);                        
                    }

                    // ... and (u\{p})
                    tempV.first = u.setMinus(p);
                    tempV.second = 0;

                    dcj.join[0] = tempU;
                    dcj.join[1] = tempV;

                    dcjs.push(dcj);
                }

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

                if(!adjA[idxU].label.empty())
                {
                    for(std::vector<int>::iterator it = adjA[idxU].label.begin();
                            it != adjA[idxU].label.end(); it++)
                    {
                        locLabelA[abs(*it)].positionLabel = idxU;
                    }
                }

                if(!adjA[idxV].label.empty())
                {
                    for(std::vector<int>::iterator it = adjA[idxV].label.begin();
                            it != adjA[idxV].label.end(); it++)
                    {
                        locLabelA[abs(*it)].positionLabel = idxV;
                    }
                }

                // Print
                Genome *g = adjacencyTableToGenome(adjA, locA);
                //printGenome(g);

                steps.push(g);

                ++dist;

            } // end if u is an adjacency
        } // end if telomere
    }// end for

    // deletion, insertion and substitution
    for(int i = 1; adjB[i].first != END_OF_TABLE; ++i)
    {
        int idxU;
        int p = adjB[i].first;

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

        if( (!adjB[i].label.empty()) && (!adjA[idxU].label.empty()) )
        {
            // Substitution
            Substitution subst;

            subst.adj = u;

            newAdjacency(idxU, i);

            // Print
            subst.print(std::cerr);

            Genome *g = adjacencyTableToGenome(adjA, locA);
            //printGenome(g);

            steps.push(g);

            dcjs.push(subst);

            ++dist;
        }

        if( (!adjB[i].label.empty()) && (adjA[idxU].label.empty()) )
        {
            // Insertion
            Insertion ins;

            ins.adj = u;

            newAdjacency(idxU, i);

            // Print

            ins.print(std::cerr);

            Genome *g = adjacencyTableToGenome(adjA, locA);
            //printGenome(g);

            steps.push(g);

            dcjs.push(ins);

            ++dist;
        }

        if( (adjB[i].label.empty()) && (!adjA[idxU].label.empty()) )
        {
            // Deletion
            Deletion del;

            del.adj = u;

            tempU.first = u.first;
            tempU.label.clear();
            tempU.second = u.second;


            // Altero a Tabela AdjA:
            adjA[idxU] = tempU;

            // Print

            del.print(std::cerr);

            Genome *g = adjacencyTableToGenome(adjA, locA);
            //printGenome(g);

            steps.push(g);

            dcjs.push(del);

            ++dist;
        }
    }
    std::cout << "Distancia: " << dist << std::endl;
    return dist;
}

int AdjacencyGraph::DCJsubstDistance(Genome *a)
{
    paths();

    int g = a->numGenes() - numLabels;
    int pL = std::min( linearSingletonInA.size(), linearSingletonInB.size() );
    int pC = std::min( circularSingletonInA.size(), circularSingletonInB.size() );
    int b = oddPaths.size() + pathsLength1.size();
    int c = cycles.size();
    int pathTable[128] = { 0 };
    int sigma = substPotential(pathTable);

    int u = getU(pathTable);
    int v = getV(pathTable);
    int w = getW(pathTable);
    int x = getX(pathTable);
    int y = getY(pathTable);
    int z = getZ(pathTable);

    int d = g - c - (b/2) + sigma - pL - pC - 2*u - 3*v - 2*w - x - 2*y - z;

    std::cout<< "g: " << g << std::endl;
    std::cout<< "c: " << c << std::endl;
    std::cout<< "b: " << b << std::endl;
    std::cout<< "sigma: " << sigma << std::endl;
    std::cout<< "pL: " << pL << std::endl;
    std::cout<< "pC: " << pC << std::endl;
    std::cout<< "u: " << u << std::endl;
    std::cout<< "v: " << v << std::endl;
    std::cout<< "w: " << w << std::endl;
    std::cout<< "x: " << x << std::endl;
    std::cout<< "y: " << y << std::endl;
    std::cout<< "z: " << z << std::endl;

    return d;
}

/*
 * Determinar o número total de runs
 */

int AdjacencyGraph::substPotential(int pathTable[128])
{
    int oddPath = 1;
    int evenPath = 2;
    int numPaths = oddPaths.size()+evenPathsFromA.size()+evenPathsFromB.size();
    
    int count = 0;

    int c = substPotentialInCycles(cycles);
    int odd = substPotentialInPaths(oddPaths, oddPath, pathTable, &count);
    int evenFromA = substPotentialInPaths(evenPathsFromA, evenPath, pathTable, &count);
    int evenFromB = substPotentialInPaths(evenPathsFromB, evenPath, pathTable, &count);

    std::cout << "potencial dos ciclos: " << c << std::endl;
    std::cout << "potencial dos caminhos ímpares: " << odd << std::endl;
    std::cout << "potencial dos caminhos pares começando em A: " << evenFromA << std::endl;
    std::cout << "potencial dos caminhos pares começando em B: " << evenFromB << std::endl;

    int substPotential = c + odd + evenFromA + evenFromB;

    for(byte i = 0; i < 128; ++i) {
        if (pathToStr(i) != NULL )
            printPath(i, std::cout) << " = " << (int)i << "\t" << pathTable[i] << std::endl;
    }

    return substPotential;

}

int AdjacencyGraph::getU(int *pathTable) {
    int u = 0;
    while (pathTable[AAab4] > 0 && pathTable[BBab4] > 0) {
        ++u;
        --pathTable[AAab4];
        --pathTable[BBab4];
    }
    return u;
}

int AdjacencyGraph::getV(int *pathTable) {
    int v = 0;
    while (pathTable[AAab4] > 1 && pathTable[BBa1] > 0 && pathTable[BBb1] > 0) {
        ++v;
        --pathTable[AAab4];
        --pathTable[AAab4];
        --pathTable[BBa1];
        --pathTable[BBb1];
    }
    while (pathTable[BBab4] > 1 && pathTable[AAa1] > 0 && pathTable[AAb1] > 0) {
        ++v;
        --pathTable[BBab4];
        --pathTable[BBab4];
        --pathTable[AAa1];
        --pathTable[AAb1];
    }
    return v;
}

int AdjacencyGraph::getW(int *pathTable) {
    int w = 0, oldw = 0;
    do {
        oldw = w;
        while (pathTable[AAab4] > 0 && pathTable[BBa1] > 0 &&
                pathTable[ABab4] > 0) {
            ++w;
            --pathTable[AAab4];
            --pathTable[BBa1];
            --pathTable[ABab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBb1] > 0 &&
                pathTable[ABba4] > 0) {
            ++w;
            --pathTable[AAab4];
            --pathTable[BBb1];
            --pathTable[ABba4];
        }
        while (pathTable[BBab4] > 0 && pathTable[AAa1] > 0 &&
                pathTable[ABba4] > 0) {
            ++w;
            --pathTable[BBab4];
            --pathTable[AAa1];
            --pathTable[ABba4];
        }
        while (pathTable[BBab4] > 0 && pathTable[AAb1] > 0 &&
                pathTable[ABab4] > 0) {
            ++w;
            --pathTable[BBab4];
            --pathTable[AAb1];
            --pathTable[ABab4];
        }
        while (pathTable[AAab4] > 1 && pathTable[BBa1] > 0) {
            ++w;
            ++pathTable[AAb3];
            --pathTable[AAab4];
            --pathTable[AAab4];
            --pathTable[BBa1];
        }
        while (pathTable[AAab4] > 1 && pathTable[BBb1] > 0) {
            ++w;
            ++pathTable[AAa3];
            --pathTable[AAab4];
            --pathTable[AAab4];
            --pathTable[BBb1];
        }
        while (pathTable[BBab4] > 1 && pathTable[AAa1] > 0) {
            ++w;
            ++pathTable[BBb3];
            --pathTable[BBab4];
            --pathTable[BBab4];
            --pathTable[AAa1];
        }
        while (pathTable[BBab4] > 1 && pathTable[AAb1] > 0) {
            ++w;
            ++pathTable[BBa3];
            --pathTable[BBab4];
            --pathTable[BBab4];
            --pathTable[AAb1];
        }
    } while (oldw != w);
    return w;
}

int AdjacencyGraph::getX(int *pathTable) {
    int x = 0, oldx = 0;
    do {
        oldx = x;
        while (pathTable[AAab4] > 1) {
            ++x;
            ++pathTable[AAa3];
            ++pathTable[AAb3];
            --pathTable[AAab4];
            --pathTable[AAab4];
        }
        while (pathTable[BBab4] > 1) {
            ++x;
            ++pathTable[BBa3];
            ++pathTable[BBb3];
            --pathTable[BBab4];
            --pathTable[BBab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[ABab4]) {
            ++x;
            ++pathTable[AAa3];
            --pathTable[AAab4];
            --pathTable[ABab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[ABba4] > 0) {
            ++x;
            ++pathTable[AAb3];
            --pathTable[AAab4];
            --pathTable[ABba4];
        }
        while (pathTable[BBab4] > 0 && pathTable[ABab4] > 0) {
            ++x;
            ++pathTable[BBb3];
            --pathTable[BBab4];
            --pathTable[ABab4];
        }
        while (pathTable[BBab4] > 0 && pathTable[ABba4] > 0) {
            ++x;
            ++pathTable[BBa3];
            --pathTable[BBab4];
            --pathTable[ABba4];
        }
        while (pathTable[AAa1] > 0 && pathTable[BBab4] > 0) {
            ++x;
            ++pathTable[ABab4];
            --pathTable[AAa1];
            --pathTable[BBab4];
        }
        while (pathTable[AAb1] > 0 && pathTable[BBab4] > 0) {
            ++x;
            ++pathTable[ABba4];
            --pathTable[AAb1];
            --pathTable[BBab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBa1] > 0) {
            ++x;
            ++pathTable[ABba4];
            --pathTable[AAab4];
            --pathTable[BBa1];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBb1] > 0) {
            ++x;
            ++pathTable[ABab4];
            --pathTable[AAab4];
            --pathTable[BBb1];
        }
        while (pathTable[AAab2] > 0 && pathTable[BBab4] > 0) {
            ++x;
            --pathTable[AAab2];
            --pathTable[BBab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBab2] > 0) {
            ++x;
            --pathTable[AAab4];
            --pathTable[BBab2];
        }
        while (pathTable[AAab2] > 0 && pathTable[BBab2] > 0) {
            ++x;
            --pathTable[AAab2];
            --pathTable[BBab2];
        }
        while (pathTable[AAa3] > 0 && pathTable[BBab4] > 0) {
            ++x;
            --pathTable[AAa3];
            --pathTable[BBab4];
        }
        while (pathTable[AAb3] > 0 && pathTable[BBab4] > 0) {
            ++x;
            --pathTable[AAb3];
            --pathTable[BBab4];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBa3] > 0) {
            ++x;
            --pathTable[AAab4];
            --pathTable[BBa3];
        }
        while (pathTable[AAab4] > 0 && pathTable[BBb3] > 0) {
            ++x;
            --pathTable[AAab4];
            --pathTable[BBb3];
        }
        while (pathTable[AAa1] > 0 && pathTable[BBa1] > 0) {
            ++x;
            --pathTable[AAa1];
            --pathTable[BBa1];
        }
        while (pathTable[AAb1] > 0 && pathTable[BBb1] > 0) {
            ++x;
            --pathTable[AAb1];
            --pathTable[BBb1];
        }
        while (pathTable[AAa1] > 0 && pathTable[BBab2] > 0) {
            ++x;
            --pathTable[AAa1];
            --pathTable[BBab2];
        }
        while (pathTable[AAb1] > 0 && pathTable[BBab2] > 0) {
            ++x;
            --pathTable[AAb1];
            --pathTable[BBab2];
        }
        while (pathTable[AAab2] > 0 && pathTable[BBa1] > 0) {
            ++x;
            --pathTable[AAab2];
            --pathTable[BBa1];
        }
        while (pathTable[AAab2] > 0 && pathTable[BBb1] > 0) {
            ++x;
            --pathTable[AAab2];
            --pathTable[BBb1];
        }
        while (pathTable[AAa1] > 0 && pathTable[BBa3] > 0) {
            ++x;
            --pathTable[AAa1];
            --pathTable[BBa3];
        }
        while (pathTable[AAb1] > 0 && pathTable[BBb3] > 0) {
            ++x;
            --pathTable[AAb1];
            --pathTable[BBb3];
        }
        while (pathTable[AAa3] > 0 && pathTable[BBa1] > 0) {
            ++x;
            --pathTable[AAa3];
            --pathTable[BBa1];
        }
        while (pathTable[AAb3] > 0 && pathTable[BBb1] > 0) {
            ++x;
            --pathTable[AAb3];
            --pathTable[BBb1];
        }
        while (pathTable[ABab4] > 0 && pathTable[ABba4] > 0) {
            ++x;
            --pathTable[ABab4];
            --pathTable[ABba4];
        }
    } while (oldx != x);
    return x;
}

int AdjacencyGraph::getY(int *pathTable) {
    int y = 0;
    while (pathTable[ABab4] > 1 && pathTable[AAb1] > 0 && pathTable[BBa1] > 0) {
        ++y;
        --pathTable[ABab4];
        --pathTable[ABab4];
        --pathTable[AAb1];
        --pathTable[BBa1];
    }
    while (pathTable[ABba4] > 1 && pathTable[AAa1] > 0 && pathTable[BBb1] > 0) {
        ++y;
        --pathTable[ABba4];
        --pathTable[ABba4];
        --pathTable[AAa1];
        --pathTable[BBb1];
    }
    return y;
}


int AdjacencyGraph::getZ(int *pathTable) {
    int z = 0, oldz = 0;
    do {
        oldz = z;
        while (pathTable[ABab4] > 0 && pathTable[AAab2] > 0 &&
                pathTable[BBa3] > 0) {
            ++z;
            --pathTable[ABab4];
            --pathTable[AAab2];
            --pathTable[BBa3];
        }
        while (pathTable[ABba4] > 0 && pathTable[AAab2] > 0 &&
                pathTable[BBb3] > 0) {
            ++z;
            --pathTable[ABba4];
            --pathTable[AAab2];
            --pathTable[BBb3];
        }
        while (pathTable[ABba4] > 0 && pathTable[AAa3] > 0 &&
                pathTable[BBab2] > 0) {
            ++z;
            --pathTable[ABba4];
            --pathTable[AAa3];
            --pathTable[BBab2];
        }
        while (pathTable[ABab4] > 0 && pathTable[AAb3] > 0 &&
                pathTable[BBab2] > 0) {
            ++z;
            --pathTable[ABab4];
            --pathTable[AAb3];
            --pathTable[BBab2];
        }
        while (pathTable[ABab4] > 0 && pathTable[AAb1] > 0 &&
                pathTable[BBa3] > 0) {
            ++z;
            --pathTable[ABab4];
            --pathTable[AAb1];
            --pathTable[BBa3];
        }
        while (pathTable[ABab4] > 0 && pathTable[AAb3] > 0 &&
                pathTable[BBa1] > 0) {
            ++z;
            --pathTable[ABab4];
            --pathTable[AAb3];
            --pathTable[BBa1];
        }
        while (pathTable[ABba4] > 0 && pathTable[AAa1] > 0 &&
                pathTable[BBb3] > 0) {
            ++z;
            --pathTable[ABba4];
            --pathTable[AAa1];
            --pathTable[BBb3];
        }
        while (pathTable[ABba4] > 0 && pathTable[AAa3] > 0 &&
                pathTable[BBb1] > 0) {
            ++z;
            --pathTable[ABba4];
            --pathTable[AAa3];
            --pathTable[BBb1];
        }
        while (pathTable[ABab4] > 0 && pathTable[AAb1] > 0 &&
                pathTable[BBa1] > 0) {
            ++z;
            ++pathTable[ABba4];
            --pathTable[ABab4];
            --pathTable[AAb1];
            --pathTable[BBa1];
        }
        while (pathTable[ABba4] > 0 && pathTable[AAa1] > 0 &&
                pathTable[BBb1] > 0) {
            ++z;
            ++pathTable[ABab4];
            --pathTable[ABba4];
            --pathTable[AAa1];
            --pathTable[BBb1];
        }
        while (pathTable[ABab4] > 1 && pathTable[AAb1] > 0) {
            ++z;
            ++pathTable[AAa3];
            --pathTable[ABab4];
            --pathTable[ABab4];
            --pathTable[AAb1];
        }
        while (pathTable[ABab4] > 1 && pathTable[BBa1] > 0) {
            ++z;
            ++pathTable[BBb3];
            --pathTable[ABab4];
            --pathTable[ABab4];
            --pathTable[BBa1];
        }
        while (pathTable[ABba4] > 1 && pathTable[AAa1] > 0) {
            ++z;
            ++pathTable[AAb3];
            --pathTable[ABba4];
            --pathTable[ABba4];
            --pathTable[AAa1];
        }
        while (pathTable[ABba4] > 1 && pathTable[BBb1] > 0) {
            ++z;
            ++pathTable[BBa3];
            --pathTable[ABba4];
            --pathTable[ABba4];
            --pathTable[BBb1];
        }
    } while (oldz != z);
    return z;
}

int AdjacencyGraph::getK(int numRuns) {
    if (numRuns == 0)
    {
    } else if (numRuns % 4 == 0) {
        return 4;
    } else {
        return numRuns % 4;
    }
}

// k = 0 means epsilon
byte AdjacencyGraph::pathToByte(WhichGenome firstElementIn,
                WhichGenome lastElementIn,
                WhichGenome firstRunIn,
                WhichGenome lastRunIn, int k) {
  byte ret = 0;

  if (firstElementIn == genomeB && lastElementIn == genomeA) {
     // replaces BA?? for AB??
    firstElementIn = genomeA;
    lastElementIn = genomeB;
    if (firstRunIn == genomeA && lastRunIn == genomeB) {
      firstRunIn = genomeB;
      lastRunIn == genomeA;
    } else if (firstRunIn == genomeB && lastRunIn == genomeA) {
      firstRunIn = genomeA;
      lastRunIn == genomeB;
    }
  } else if (firstElementIn == lastElementIn) {
    // replaces AAba+2, AAba+4, BBba+2, BBba+4 for
    //          AAab+2, AAab+4, BBab+2, BBab+4
    if (firstRunIn == genomeB && lastRunIn == genomeA) {
      firstRunIn = genomeA;
      lastRunIn = genomeB;
    }
  }

  if (firstElementIn == genomeB) {
    ret |= 64; // 64 = 1000000b
  }
  if (lastElementIn == genomeB) {
    ret |= 32; // 32 = 100000b
  }
  if (k > 0) {
    if (firstRunIn == genomeB) {
      ret |= 16; // 16 = 10000b
    }
    if (k%2 == 0) {
      if (lastRunIn == genomeB) {
        ret |= 8; // 16 = 1000b
      }
    }
    ret |= k & 7 ; //  K & 111b == k2 k1 k0
  }
  return ret;
}

const char *AdjacencyGraph::pathToStr(byte representation) {
    static const
    char *strByte[] = { "AAe"   , "AAa+1" , NULL     , "AAa+3" , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 10 */   "AAab+2", NULL    , "AAab+4" , NULL    , NULL     ,
                        NULL    , NULL    , "AAb+1"  , NULL    , "AAb+3"  ,
             /* 20 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 30 */   NULL    , NULL    , "ABe"    , "ABa+1" , NULL     ,
                        "ABa+3" , NULL    , NULL     , NULL    , NULL     ,
             /* 40 */   NULL    , NULL    , "ABab+2" , NULL    , "ABab+4" ,
                        NULL    , NULL    , NULL     , NULL    , "ABb+1"  ,
             /* 50 */   "ABba+2", "ABb+3" , "ABba+4" , NULL    , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 60 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 70 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 80 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , NULL    , NULL     , NULL    , NULL     ,
             /* 90 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , "BBe"   , "BBa+1"  , NULL    , "BBa+3"  ,
             /*100 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , "BBab+2", NULL     , "BBab+4", NULL     ,
             /*110 */   NULL    , NULL    , NULL     , "BBb+1" , NULL     ,
                        "BBb+3" , NULL    , NULL     , NULL    , NULL     ,
             /*120 */   NULL    , NULL    , NULL     , NULL    , NULL     ,
                        NULL    , NULL    , NULL };
    return strByte[representation];
}

std::ostream &AdjacencyGraph::printPath(byte representation,
                                        std::ostream &os) {
      const char *str = pathToStr(representation);
      if (str != NULL)
          os << str;
      else
          os << "invalid";
      return os;
  }

/**
 * Encontra o índice de uma extremidade de gene na tabela de adjacência,
 * dados o gene e a tabela de localização correspondente.
 */
int find(int geneEnd, Location *locTable) {
        if(geneEnd > 0)
            return locTable[geneEnd].tail;
        else if(geneEnd < 0)
            return locTable[-geneEnd].head;
}

/*
 * Determinar o número de runs nos ciclos
 */
int AdjacencyGraph::substPotentialInCycles(std::deque< std::pair<WhichGenome,int> > cycle)
{
    int substPTotal = 0; // Somatório de potenciais de substituição

    std::deque< std::pair<WhichGenome,int> >::iterator it;

    for(it = cycle.begin(); it != cycle.end(); it++)
    {
        int numRuns = 0;
        std::pair<WhichGenome,int> p = *it;

        WhichGenome firstGenome = p.first, genome = p.first;
        int firstGene = p.second, gene = p.second;

        WhichGenome lastLabelIn = undef;

        do {
          Adjacency *adjTable;

          int idx;
          if (genome == genomeA) {
            adjTable = adjA;
            idx = find(gene, locA);
          } else {
            adjTable = adjB;
            idx = find(gene, locB);
          }

          if (!adjTable[idx].label.empty() && lastLabelIn != genome) {
            ++numRuns;
            lastLabelIn = genome;
          }
          if (genome == genomeA) {
            genome = genomeB;
          } else {
            genome = genomeA;
          }
          gene = adjTable[idx].setMinus(gene);
        } while(gene != firstGene || genome != firstGenome);

        if (numRuns > 1 && ((numRuns % 2) == 1)) {
          --numRuns;
        }

        int substP; // potencial de substituição da componente c
        if(numRuns >= 1)
            substP = ceil((numRuns+1.0)/4.0);
        else
            substP = 0;

        substPTotal = substPTotal + substP;
    }
    return substPTotal;
}

/*
 * Determinar o número de runs nos caminhos pares e ímpares
 */

int AdjacencyGraph::substPotentialInPaths(std::deque<Path> paths,
                                int parity, int *pathTable, int *count)
{
    WhichGenome lastLabelIn = undef;
    WhichGenome whereThis = undef;

    WhichGenome startsRunIn = undef;
    WhichGenome endsRunIn = undef;

    WhichGenome startsPathIn = undef;
    WhichGenome endsPathIn = undef;

    int k;

    int aRuns, bRuns, numRuns;
    int i, idx, temp;
    int first; // Condição de parada do ciclo
    
    Adjacency *adjTable;
    Location *locTable;
    
    int substPTotal = 0; // Somatório de potenciais de substituição
    int substP; // potencial de substituição da componente c

    std::deque<Path>::iterator it;

    // Determinar o número de runs

    for(it = paths.begin(); it != paths.end(); it++)
    {
        aRuns = 0;
        bRuns = 0;
        numRuns = 0;
        Path p = *it;

        i = p.start;

        if(parity == 1)
        {
            // Odd Path
            if(p.startsIn == genomeA)
            {
                startsPathIn = genomeA;
                endsPathIn = genomeB;
            }
            else if(p.startsIn == genomeB)
            {
                startsPathIn = genomeB;
                endsPathIn = genomeA;
            }
        }
        else if(parity == 2)
        {
            // Even Path
            if(p.startsIn == genomeA)
            {
                startsPathIn = genomeA;
                endsPathIn = genomeA;
            }
            else if(p.startsIn == genomeB)
            {
                startsPathIn = genomeB;
                endsPathIn = genomeB;
            }
        }

        if(p.startsIn == genomeA)
        {
            whereThis = genomeA;
            adjTable = adjA;
            locTable = locA;
        }
        else if(p.startsIn == genomeB)
        {
            whereThis = genomeB;
            adjTable = adjB;
            locTable = locB;
        }

        if(i > 0)
            idx = locTable[i].tail;
        else if(i < 0)
            idx = locTable[-i].head;

        for(;;)
        {   
            // A adjacencia possui label?
            if(!adjTable[idx].label.empty())
            {
                if(startsRunIn == undef)
                    startsRunIn = whereThis;

                if(lastLabelIn != whereThis)
                {
                    if(whereThis == genomeA)
                        ++aRuns;
                    if(whereThis == genomeB)
                        ++bRuns;
                }
                lastLabelIn = whereThis;
            }

            // Condição de parada
            if(i == 0)
            {
                endsRunIn = lastLabelIn;
                break;
            }
            // Atualizando a adjacencia
            if(whereThis == genomeA)
            {
                whereThis = genomeB;
                adjTable = adjB;
                locTable = locB;
            }
            else if(whereThis == genomeB)
            {
                whereThis = genomeA;
                adjTable = adjA;
                locTable = locA;
            }

            if(i > 0)
                idx = locTable[i].tail;
            else if(i < 0)
                idx = locTable[-i].head;

            temp = i;

            i = adjTable[idx].setMinus(temp);
        }

        numRuns = aRuns + bRuns;

        if(numRuns >= 1)
            substP = ceil((numRuns+1.0)/4.0);
        else
            substP = 0;

        substPTotal = substPTotal + substP;

        k = getK(numRuns);

        byte type = pathToByte(startsPathIn, endsPathIn, startsRunIn, endsRunIn, k);
        pathTable[type]++;
    }
    return substPTotal;
}

void AdjacencyGraph::findLabels(Genome *a, Genome *b)
{
    labels.clear();
    numLabels = 0;

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
            {
                labelsInA[k].insert(abs((*chrA)[i]));
                numLabels++;
            }
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
            {
                labelsInB[k].insert(abs((*chrB)[i]));
            }
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
*/

std::ostream & operator<<(std::ostream &os, const std::vector<int> &l) {
    std::vector<int>::const_iterator _it = l.begin();
    os << "<";
    while (_it != l.end()) {
        if (_it != l.begin()) {
            os << ",";
        }
        os << *_it;
        ++_it;
    }
    os << ">";
    return os;
}

