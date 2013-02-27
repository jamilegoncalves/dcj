#include "AdjacencyGraph.h"

AdjacencyGraph::AdjacencyGraph(Genome *a, Genome *b)
{
    // Bergeron, Mixtacki, Stoye. A Unifying View
    // of Genome earrangements. 2006. Algorithm 1
    n = a->numGenes();

    vertex = new Vertex[2*n + 1];

    constructTables(a, adjA, locA);
    constructTables(b, adjB, locB);

    adjA[10].first = 0;
    adjB[9].first = 0;
    
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
    
    for(int i = 1; i <= 15; ++i)
    {
        std::cerr << "AdjA["<<i<<"].first: "<< adjA[i].first << std::endl;
        std::cerr << "AdjA["<<i<<"].second: "<< adjA[i].second << std::endl;
    }
    
    for(int i = 1; i <= 14; ++i)
    {
        std::cerr << "AdjB["<<i<<"].first: "<< adjB[i].first << std::endl;
        std::cerr << "AdjB["<<i<<"].second: "<< adjB[i].second << std::endl;
    }
}

AdjacencyGraph::~AdjacencyGraph()
{
    delete vertex;
    delete adjA, adjB, locA, locB;
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
            adj[offset + 1].first = chr->genes[1];
            adj[offset + 1].second = 0;
            adj[offset + chr->genes.size() +1].first = - chr->genes[chr->genes.size()];
            adj[offset + chr->genes.size() +1].second = 0;
            
            for(int i = 2; i <= chr->genes.size(); ++i)
            {
                adj[offset + i].first = - chr->genes[i-1];
                adj[offset + i].second = chr->genes[i];
            }
        }
        else // if it's a circular chromosome
        {
            if(chr->genes.size() == 1)
            {
                adj[offset + 1].first = chr->genes[1];
                adj[offset + 1].second = chr->genes[1];
            }
            else
            {
                adj[offset + 1].first = chr->genes[0];
                adj[offset + chr->genes.size()].second = - chr->genes[0];
                
                for(int i = 1; i < chr->genes.size(); ++i)
                {
                    adj[offset + i].second = - chr->genes[i];
                    adj[offset + i + 1].first = chr->genes[i];
                }
            }
        }
        
        // Constroi tabela 2:
        for(int i = 1; i <= chr->genes.size(); ++i)
        {
            int label = adj[offset + i].first;
        
            if(label > 0)
                loc[label].tail = i;
            else if(label < 0)
                loc[-label].head = i;
        
            label = adj[offset + i].second;
        
            if(label == 0)
                continue;
            else if(label > 0)
                loc[label].tail = label;
            else if(label < 0)
                loc[-label].head = -label;
        }
        offset = offset + chr->numAdjacencies(chr);
    }
    
/*
    for(int i = 0; i <= offset; ++i)
    {
        std::cout << "Adj["<< i <<"].first: " << adj[i].first << std::endl;
        std::cout << "Adj["<< i <<"].second: " << adj[i].second << std::endl;
        
        std::cout << "Loc["<< i <<"].head: " << loc[i].head << std::endl;
        std::cout << "Loc["<< i <<"].tail: " << loc[i].tail << std::endl;        
    }
 */
}

int AdjacencyGraph::numCycles()
{
    
}

int AdjacencyGraph::oddPaths()
{
    
}

int AdjacencyGraph::dcjDistance(Genome *a, Genome *b)
{
    int n = a->numGenes();
    int c = numCycles();
    int i = oddPaths();
    
    int dist = n - (c + i/2);

    return dist;
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
                newadjA[u] = adjB[i];
                newadjA[v].first = adjA[u].second;
                newadjA[v].second = adjA[v].second;

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
        
        int n = chr->genes.size();
        
        if(chr->isLinear() == true)
            totalAdj = totalAdj + n+1;
        else
            totalAdj = totalAdj + n;
    }
    return totalAdj;
}


std::ostream & operator<<(std::ostream &os, const AdjacencyGraph& ag)
{
    for (int i = 1; i <= ag.n; ++i)
        os << ag.vertex[i].vertex1 << "/"
           << ag.vertex[i].vertex2 << std::endl;

    return os;
}