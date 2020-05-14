#include <iostream>

#include "Grid.h"
#include "UniversalElement.h"

int main() {

    Grid grid;
    grid.printGrid();

    UniversalElement uelem;


    uelem.createMatrixHandC(grid.elements[0]);
    uelem.print();

    grid.aggregationHandC();
    grid.egdeLength();

    //uelem.vectorP();

}