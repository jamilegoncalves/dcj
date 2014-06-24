#include "Rearrangement.h"

void Rearrangement::print(std::ostream &os) const
{
}

std::ostream & operator<<(std::ostream &os, const Rearrangement &dcj)
{
    dcj.print(os);
    return os;
}

