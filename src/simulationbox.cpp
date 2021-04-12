#include "simulationbox.h"
#include "constants.h"

Box::Box(double inx, double iny, double inz){
    /* Construct the simulation box, converting to
     * atomic units as we go */
    sides << inx * constants::bohrAng, iny * constants::bohrAng, inz * constants::bohrAng;
    dimensions = 3;
}

double Box::volume(){
    return sides(0) * sides(1) * sides(2);
}