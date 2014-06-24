/* 
 * File:   Insertion.h
 * Author: jamile
 *
 * Created on 28 de Abril de 2014, 18:22
 */

#ifndef INSERTION_H
#define	INSERTION_H

#include "Adjacency.h"
#include "Rearrangement.h"
#include <queue>

class Insertion : public Rearrangement {
public:
    Insertion();
    void print(std::ostream &os) const;

    Adjacency adj;
    std::vector<int> label;
};

#endif	/* INSERTION_H */

