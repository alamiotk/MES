//
// Created by ala on 06.05.20.
//
#include "UniversalElement.h"
#include <math.h>
#include <iostream>
#include "Element.h"


using namespace std;

const double VALUE_PLUS = (1/(sqrt(3)));
const double VALUE_MINUS = - VALUE_PLUS;

UniversalElement::UniversalElement() {
    this -> integrationPoints[0][0] = VALUE_MINUS;
    this -> integrationPoints[0][1] = VALUE_MINUS;
    this -> integrationPoints[1][0] = VALUE_PLUS;
    this -> integrationPoints[1][1] = VALUE_MINUS;
    this -> integrationPoints[2][0] = VALUE_PLUS;
    this -> integrationPoints[2][1] = VALUE_PLUS;
    this -> integrationPoints[3][0] = VALUE_MINUS;
    this -> integrationPoints[3][1] = VALUE_PLUS;

    this -> calculateShapeFunctions();

}

UniversalElement::~UniversalElement() {
    // delete[] nd;

}


void UniversalElement::calculateShapeFunctions() {
    for (int i = 0; i < globalData.numberOfNodesInElement; i++) {
        for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
            shapeFunctionsKsiEta(this->integrationPoints[i][0], this->integrationPoints[i][1]);
            derativeKsi[i][j] = dNdKsi[j];                                                      // derivates ksi
            derativeEta[i][j] = dNdEta[j];                                                      // derivates eta

            shapeFunctionsN(this -> integrationPoints[i][0], this -> integrationPoints[i][1]);
            functionsN[i][j] = N[j];                                                            // functions N
        }
    }
}


void UniversalElement::shapeFunctionsN(double ksi, double eta) {
    N[0] = 0.25 * (1 - ksi) * (1 - eta);
    N[1] = 0.25 * (1 + ksi) * (1 - eta);
    N[2] = 0.25 * (1 + ksi) * (1 + eta);
    N[3] = 0.25 * (1 - ksi) * (1 + eta);
}



void UniversalElement::shapeFunctionsKsiEta(double ksi, double eta) {
    dNdKsi[0] = -(0.25 * (1 - eta));
    dNdEta[0] = -(0.25 * (1 - ksi));

    dNdKsi[1] = (0.25 * (1 - eta));
    dNdEta[1] = -(0.25*(1 + ksi));

    dNdKsi[2] = (0.25 * (1 + eta));
    dNdEta[2] = (0.25 * (1 + ksi));

    dNdKsi[3] = -(0.25 * (1 + eta));
    dNdEta[3] = (0.25 * (1 - ksi));
}


void UniversalElement::calculcateJacobiTransformation(Element element) {
//    double x[] = {0, 0.025, 0.025, 0};
//    double y[] = {0, 0, 0.025, 0.025};
//
//    for ( int i = 0; i < globalData.numberOfNodesInElement; i++) {
//        jacobiMatrix[0][0] += derativeKsi[0][i] * x[i];
//        jacobiMatrix[0][1] += derativeEta[0][i] * x[i];
//
//        jacobiMatrix[1][0] += derativeKsi[0][i] * y[i];
//        jacobiMatrix[1][1] += derativeEta[0][i] * y[i];
//
//    }

    for ( int i = 0; i < globalData.numberOfNodesInElement; i++) {
        jacobiMatrix[0][0] += derativeKsi[0][i] * element.nodes[i]->x;
        jacobiMatrix[0][1] += derativeEta[0][i] * element.nodes[i]->x;

        jacobiMatrix[1][0] += derativeKsi[0][i] * element.nodes[i]->y;
        jacobiMatrix[1][1] += derativeEta[0][i] * element.nodes[i]->y;

    }

    det = jacobiMatrix[0][0]*jacobiMatrix[1][1] - jacobiMatrix[0][1]*jacobiMatrix[1][0];

    cout << det<< endl;

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            jacobiReverseMatrix[i][j] = jacobiMatrix[i][j]*(1/det);     // jacobi's reverse matrix
        }
    }
}

void UniversalElement::createMatrixHandC(Element element) {

    this -> calculcateJacobiTransformation(element);

    for (int i = 0; i < globalData.numberOfNodesInElement; i++) {
        for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
            derivateX[i][j] =(
                jacobiReverseMatrix[0][0] * derativeKsi[i][j] +
                jacobiReverseMatrix[0][1] * derativeEta[i][j]
            );
            derivateY[i][j] = (
                jacobiReverseMatrix[1][0] * derativeKsi[i][j] +
                jacobiReverseMatrix[1][1] * derativeEta[i][j]
            );
        }
    }

    for (int i = 0 ; i < 4; i++) {
        for (int j=0; j< 4; j++ ){
            H[i][j] = 0;
            C[i][j] = 0;
        }
    }
//    cout << "Stałe" << globalData.heat << " " << globalData.density << " " << det << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            for ( int m = 0; m < 4; m++){
                H[j][m] += (derivateX[i][j]*derivateX[i][m] + derivateY[i][j]*derivateY[i][m])*globalData.conductivity*det;

//                cout << functionsN[i][j] << " " << functionsN[i][m] << " " << (functionsN[i][j]*functionsN[i][m]) << endl;
//                cout << C[j][m] << endl;
                C[j][m] += (functionsN[i][j]*functionsN[i][m])*globalData.heat*globalData.density*det;
            }
        }
    }
}

void UniversalElement::print() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {                                   // derivated ksi
            cout << derativeKsi[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {                                   // derivates eta
            cout << derativeEta[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {                                   // jacobi's matrix
            cout << jacobiMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {                                   // jacobi's reverse matrix
            cout << jacobiReverseMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {                           // matrix H
            cout << H[i][j] << "  ";
        }
        cout << endl;
    }


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {                                   // functionsN
            cout << functionsN[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {                           // matrix C
            cout << C[i][j] << "  ";
        }
        cout << endl;
    }

}
