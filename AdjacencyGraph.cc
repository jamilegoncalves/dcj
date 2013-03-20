#include "AdjacencyGraph.h"

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Bergeron, Mixtacki, Stoye. A Unifying View
    // of Genome earrangements. 2006. Algorithm 1
    n = a->numGenes();
    numAdj = totalAdjacencies(a);
    vertex = new Vertex[2*n + 1];

    adjA = adjB = NULL; // DEBUG

    constructTables(a, adjA, locA);
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

void AdjacencyGraph::printAdjacencies(std::ostream &os)
{
    if ( adjA != NULL )
    {
        os << "AdjA:" << std::endl;
        os << "first:";
        for(int i = 1; i <= 10; ++i)
            os << "\t" << adjA[i].first;
        os << std::endl;

        os << "second:";
        for(int i = 1; i <= 9; ++i)
            os << "\t" << adjA[i].second;
        os << std::endl;
    }

    if ( locA != NULL )
    {
        os << "LocA:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 7; ++i)
            os << "\t" << locA[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 7; ++i)
            os << "\t" << locA[i].tail;
        os << std::endl;
    }

    if ( adjB != NULL )
    {
        os << std::endl;
        os << "AdjB:" << std::endl;
        os << "first:";
        for(int i = 1; i <= 9; ++i)
            os << "\t" << adjB[i].first;
        os << std::endl;

        os << "second:";
        for(int i = 1; i <= 9; ++i)
            os << "\t" << adjB[i].second;
        os << std::endl;
    }

    if ( locB != NULL )
    {
        os << std::endl;
        os << "LocB:" << std::endl;
        os << "head:";
        for(int i = 1; i <= 7; ++i)
            os << "\t" << locB[i].head;
        os << std::endl;

        os << "tail:";
        for(int i = 1; i <= 7; ++i)
            os << "\t" << locB[i].tail;
        os << std::endl;
    }
}

int AdjacencyGraph::sortByDCJ()
{
    return sortByDCJ(n, adjA, adjB, locA, locB);
}

/**
 * Constroi tabelas: Adjacency e Location
 */

void AdjacencyGraph::constructTables(Genome *g, Adjacency *&adj, Location *&loc)
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
                adj[offset + chr->length()].second = - (*chr)[1];

                for(int i = 2; i <= chr->length(); ++i)
                {
                    adj[offset + 1].second = - (*chr)[i];
                    adj[offset + 2].first = (*chr)[i];
                }
            }
        }

        // Constroi tabela 2:
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
                loc[-label].head = i;
        }
        offset = offset + chr->numAdjacencies(chr);
    } // end for
    adj[totalAdjacencies(g)+1].first = 0;
}

/**
 * Implementation of Bergeron, Mixtacki & Stoye's greedy sorting
 * by DCJ (Algorithm 2).
 */
int AdjacencyGraph::sortByDCJ(int n, Adjacency *adjA, Adjacency *adjB,
                                    Location *locA, Location *locB)
{
    Adjacency *newadjA = new Adjacency[2*n + 2];

    int dist = 0;

    // iterate over all adjacencies/telomeres of genome B
    for(int i = 1; adjB[i].first != 0; ++i)
    {
        int u, v;

        if(adjB[i].second != 0) // if it is an adjacency
        {
            if(adjB[i].first > 0)
                u = locA[adjB[i].first].tail;
            else
                u = locA[- adjB[i].first].head;

            if(adjB[i].second > 0)
                v = locA[adjB[i].second].tail;
            else
                v = locA[- adjB[i].second].head;

            if(u != v)
            {
                std::cout << "Cut: " << adjA[u].first << "," << adjA[u].second
                                                << std::endl;
                std::cout << "Cut: " << adjA[v].first << "," << adjA[v].second
                                                << std::endl;

                newadjA[u] = adjB[i];
                newadjA[v].first = adjA[u].second;
                newadjA[v].second = adjA[v].second;

                std::cout << "Join: " << newadjA[u].first << ","
                                            << newadjA[u].second << std::endl;
                std::cout << "Join: " << newadjA[v].first << ","
                                            << newadjA[v].second << std::endl;

                // Altero a Tabela AdjA:
                adjA[u] = newadjA[u];
                adjA[v] = newadjA[v];

                // Altero a Tabela LocA:
                if(adjA[u].first > 0)
                    locA[adjA[u].first].tail = u;
                else
                    locA[-adjA[u].first].head = u;

                if(adjA[u].second > 0)
                    locA[adjA[u].second].tail = u;
                else
                    locA[-adjA[u].second].head = u;

                if(adjA[v].first > 0)
                    locA[adjA[v].first].tail = v;
                else
                    locA[-adjA[v].first].head = v;

                if(adjA[v].second > 0)
                    locA[adjA[v].second].tail = v;
                else
                    locA[-adjA[v].second].head = v;

                ++dist;
            }
        }
        else // if it is a telomere in B
        {
            if(adjB[i].first > 0)
                u = locA[adjB[i].first].tail;
            else
                u = locA[- adjB[i].first].head;

            if(adjA[u].second != 0) // if u is an adjacency
            {
                newadjA[u].first = adjB[i].first;
                newadjA[u].second = adjB[u].second;

                // Altero Tabela AdjA:
                adjA[u] = newadjA[u];

                // Altero Tabela LocA:
                if(adjA[u].first > 0)
                    locA[adjA[u].first].tail = u;
                else
                    locA[-adjA[u].first].head = u;

                if(adjA[u].second > 0)
                    locA[adjA[u].second].tail = u;
                else
                    locA[-adjA[u].second].head = u;

                ++dist;
            } // end if u is an adjacency

        } // end if adjacency/telomere
    } // end for

    delete newadjA;
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


std::ostream & operator<<(std::ostream &os, const AdjacencyGraph& ag)
{
    for (int i = 1; i <= ag.numAdj; ++i)
        os << ag.vertex[i].vertex1 << "/"
           << ag.vertex[i].vertex2 << std::endl;

    return os;
}