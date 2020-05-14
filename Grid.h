//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_GRID_H
#define MES_PROJ_GRID_H

#include "vector"

#include "GlobalData.h"
#include "Element.h"
#include "Node.h"
#include "UniversalElement.h"

using namespace std;


class Grid: public GlobalData {
public:
    vector<Element> elements;
    vector<Node> nodes;

    vector<vector<double> > aggregationMatrixH;
    vector<vector<double> > aggregationMatrixC;



    double deltaX, deltaY;

    Grid ();
    bool checkBorderCondition(double x, double y);
    void printGrid();
    void aggregationHandC();
    void egdeLength();

};


#endif //MES_PROJ_GRID_H
