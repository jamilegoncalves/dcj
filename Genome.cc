/* Genome.cc - Class for representing a genome
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
#include "Genome.h"

Genome::Genome(std::string name)
{
    this->name = name;
}

/**
 * Determinar número de genes no genoma
 */

int Genome::numGenes()
{
    int numGenes = 0;

    std::vector<Chromosome*>::iterator cIterator;

    for(cIterator = chromosomes.begin();
            cIterator != chromosomes.end(); ++cIterator)
    {
        Chromosome *chr = *cIterator;

        numGenes = numGenes + chr->length();
    }
    return numGenes;
}

std::ostream & operator<<(std::ostream &os, const Genome& g)
{
    os << ">" << g.name << std::endl;

    for ( int i = 0; i < g.chromosomes.size(); ++i )
        os << *(g.chromosomes[i]);
    os << std::endl;

    return os;
}