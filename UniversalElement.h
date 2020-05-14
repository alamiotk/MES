//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_UNIVERSALELEMENT_H
#define MES_PROJ_UNIVERSALELEMENT_H

#include "Element.h"
#include "GlobalData.h"

class UniversalElement {
public:
    GlobalData globalData;

    double integrationPoints[4][2];
    double pointsEdges[8][2];
    double dNdKsi[4];
    double dNdEta[4];
    double N[4];
    double jacobiReverseMatrix[2][2];
    double jacobiMatrix[2][2];
    double det;
    double H[4][4] = {};
    double C[4][4] = {};
    double HB[4][4] = {};
    double functionsN[4][4];
    double derivativeDNDKsi[4][4];
    double derivativeDNDEta[4][4];
    double derivativeDNDX[4][4];
    double derivativeDNDY[4][4];

    double functionNVectoP[8][4] = {};

//    Element *element;

    UniversalElement();
    ~UniversalElement();
    void calculateShapeFunctions();
    void shapeFunctionsN(double ksi, double eta);
    void shapeFunctionsKsiEta(double ksi, double eta);
    void calculcateJacobiTransformation(Element element);
    void createMatrixHandC(Element element);
    void pointsOnTheEdges();
    void print();


    void vectorP(double detJ);

};


#endif //MES_PROJ_UNIVERSALELEMENT_H