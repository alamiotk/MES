//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_UNIVERSALELEMENT_H
#define MES_PROJ_UNIVERSALELEMENT_H

#include "Element.h"
#include "GlobalData.h"
#include <vector>

using namespace std;


class UniversalElement {
public:
    GlobalData globalData;

    double jacobiReverseMatrix[2][2];
    double HBC[4][4] = {};
    double vecP[4] = {};
    double jacobiMatrix[2][2];
    double det;

    vector<vector<double> > integrationPoints;
    vector<vector<double> > pointsEdges;
    vector<double> dNdKsi;
    vector<double> dNdEta;
    vector<double> N;
    vector<vector<double> > H;
    vector<vector<double> > C;
    vector<vector<double> > derivativeDNDX;
    vector<vector<double> > derivativeDNDY;
    vector<vector<double> > functionsN;
    vector<vector<double> > functionNEdges;
    vector<vector<double> > derivativeDNDKsi;
    vector<vector<double> > derivativeDNDEta;




    UniversalElement();
    void calculateShapeFunctions();
    void shapeFunctionsN(double ksi, double eta);
    void shapeFunctionsKsiEta(double ksi, double eta);
    void calculcateJacobiTransformation(Element element, int p);
    void createMatrixHandC(Element element);
    void pointsOnTheEdges();
    void printLocal();



    void matrixHBCandVecP(Element elements, Node *node1, Node *node2, double detJ);

};


#endif //MES_PROJ_UNIVERSALELEMENT_H
