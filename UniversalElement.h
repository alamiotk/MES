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
    double dNdKsi[4];
    double dNdEta[4];
    double N[4];
    double jacobiReverseMatrix[4][4];
    double jacobiMatrix[2][2] = {};
    double det;
    double H[4][4] = {};
    double C[4][4] = {};
    double functionsN[4][4];
    double derativeKsi[4][4];
    double derativeEta[4][4];
    double derivateX[4][4];
    double derivateY[4][4];

//    Element *element;

    UniversalElement();
    ~UniversalElement();
    void calculateShapeFunctions();
    void shapeFunctionsN(double ksi, double eta);
    void shapeFunctionsKsiEta(double ksi, double eta);
    void calculcateJacobiTransformation(Element element);
    void createMatrixHandC(Element element);
    void print();

};


#endif //MES_PROJ_UNIVERSALELEMENT_H
