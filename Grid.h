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
    vector<vector<double> > aggregationMatrixHBC;



    double deltaX, deltaY;
    double detJ;

    Grid ();
    bool checkBorderCondition(double x, double y);
    void printGrid();
    void aggregationHandC();
    double edgeLength(Node *nodes1, Node *nodes2);

    bool checkIfEdge(Element elements);

};


#endif //MES_PROJ_GRID_H
