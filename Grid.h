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

    double deltaX, deltaY;
    double detJ;
    unsigned int aggregationMatrixSize;
    vector<vector<double> > aggregationMatrixH;
    vector<vector<double> > aggregationMatrixC;
    vector<double> aggregationVectoP;
    vector<double> temperatureInitialMatrix;
    vector<double>  matrixCxT0;
    vector<double> temperatureT1;


    Grid ();
    bool checkBorderCondition(double x, double y);
    void aggregation();
    void checkIfEdge(Element elements, UniversalElement universalElement, double vectorPLocalResult[],double matrixHBLocalResult[][4]);
    double edgeLength(Element elements, Node *nodes1, Node *nodes2);
    double min(vector<double> temperatureT1, int aggregationMatrixSize);
    double max(vector<double> temperatureT1, int aggregationMatrixSize);
    vector<double> solveEquation(vector<vector<double> >aggregationMatrixH, vector<double> aggregationVectoP, int aggregationMatrixSize);

    void printGrid();
    void aggregationPrint();
    void printVectorP();
    void printTemperatureT1(double minimalTemperature, double maximalTemperature);
};


#endif //MES_PROJ_GRID_H
