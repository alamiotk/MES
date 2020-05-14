#include <iostream>

#include "Grid.h"
#include "UniversalElement.h"

int main() {

    Grid grid;
    grid.printGrid();

    UniversalElement uelem;


    uelem.createMatrixHandC(grid.elements[0]);

    grid.aggregationHandC();
    uelem.print();

    //grid.egdeLength();

    //uelem.vectorP();

}