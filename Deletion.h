/* 
 * File:   Deletion.h
 * Author: jamile
 *
 * Created on 28 de Abril de 2014, 18:22
 */

#ifndef DELETION_H
#define	DELETION_H

#include "Adjacency.h"
#include "Rearrangement.h"
#include <queue>

class Deletion : public Rearrangement {
public:
    Deletion();
    void print(std::ostream &os) const;

    Adjacency adj;
};

#endif	/* DELETION_H */

