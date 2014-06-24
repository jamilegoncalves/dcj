#include "DoubleCutAndJoin.h"
#include "AdjacencyGraph.h"
#include <iostream>
#include <queue>

DoubleCutAndJoin::DoubleCutAndJoin()
{

}

void DoubleCutAndJoin::print(std::ostream &os) const
{
    os << "\n dcj cut: (" << cut[0].first << "," << cut[0].second <<
          ") and (" << cut[1].first << "," << cut[1].second << ")";

    os << "dcj join: (" << join[0].first << "," << join[0].second <<
          ") and (" << join[1].first << "," << join[1].second << ")";
    os << "\n";
    os << "\n";
}