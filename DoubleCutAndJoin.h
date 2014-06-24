/* 
 * File:   DoubleCutAndJoin.h
 * Author: jamile
 *
 * Created on 14 de Abril de 2014, 17:14
 */

#ifndef DOUBLECUTANDJOIN_H
#define	DOUBLECUTANDJOIN_H

#include "Adjacency.h"
#include "Genome.h"
#include "Rearrangement.h"
#include <queue>

class DoubleCutAndJoin : public Rearrangement {
public:
    DoubleCutAndJoin();
    void print(std::ostream &os) const;

    Adjacency cut[2];
    Adjacency join[2];

};

#endif	/* DOUBLECUTANDJOIN_H */

