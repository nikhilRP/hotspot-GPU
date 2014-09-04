#ifndef HOTSPOT_HPP_
#define HOTSPOT_HPP_

#include <cstdio>

namespace nvbio {
namespace hotspot {

struct Hotspot
{
    int filterDist;
    int filterIndexLeft;
    int filterIndexRight;
    int filterDensIndexLeft;
    int filterDensIndexRight;
    int filterSize;
    double filterWidth;
    int minSite;
    int maxSite;
    int averagePos;
    int maxWindow;

    // Statistics
    double filteredZScore;
    double filteredZScoreAdjusted;
    double weightedAvgSD;

    // Background data
    double densCount;

    /**
    * Initialize a Hotspot
    */
    Hotspot()
        : filterDist( -1 ),
        filterIndexLeft( -1 ),
        filterIndexRight( -1 ),
        filterDensIndexLeft( -1 ),
        filterDensIndexRight( -1 ),
        filterSize( 0 ),
        filterWidth( 0.0 ),
        minSite( -1 ),
        maxSite( -1 ),
        averagePos( -1 ),
        maxWindow( -1 ),
        filteredZScore( 0.0 ),
        filteredZScoreAdjusted( 0.0 ),
        weightedAvgSD( 0.0 ),
        densCount( 0.0 )
    { /* */ }
};

} // hotspot
} // nvbio
#endif /* HOTSPOT_HPP_ */
