/* 
 * File:   Rearrangement.h
 * Author: jamile
 *
 * Created on 28 de Abril de 2014, 18:12
 */

#ifndef REARRANGEMENT_H
#define	REARRANGEMENT_H
#include <iostream>

class Rearrangement
{
    friend std::ostream & operator<<(std::ostream &os,
                                     const Rearrangement &dcj);
    virtual void print(std::ostream &os) const;
};

std::ostream & operator<<(std::ostream &os, const Rearrangement &dcj);
#endif	/* REARRANGEMENT_H */

