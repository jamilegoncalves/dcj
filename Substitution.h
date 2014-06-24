/* 
 * File:   Substitution.h
 * Author: jamile
 *
 * Created on 28 de Abril de 2014, 18:23
 */

#ifndef SUBSTITUTION_H
#define	SUBSTITUTION_H

#include "Adjacency.h"
#include "Rearrangement.h"
#include <queue>

class Substitution : public Rearrangement {
public:
    Substitution();
    void print(std::ostream &os) const;

    Adjacency adj;
    std::vector<int> label;
};

#endif	/* SUBSTITUTION_H */

