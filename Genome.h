/* Genome.h - Class for representing a genome
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
#ifndef GENOME_H
#define GENOME_H
#include "Chromosome.h"
#include <vector>
#include <string>
#include <iostream>

class Genome
{
public:
	//Genome(std::string name);
        Genome(std::string name, std::string description);
        int numGenes();

        /**
         * Prints the representation of this Genome to a
         * stream.
         */
        friend std::ostream & operator<<( std::ostream &os,
                                          const Genome& g);

	std::vector<Chromosome *> chromosomes;
	std::string name;
};

#endif

