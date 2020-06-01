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

//------------------INTEGRATION POINTS ON 2D IN UNIVERSAL ELEMENT------------------------------
UniversalElement::UniversalElement() {
    this -> dNdKsi.resize(
        globalData.numberOfNodesInElement, 0
    );
    this -> dNdEta.resize(
        globalData.numberOfNodesInElement, 0
    );
    this -> N.resize(
        globalData.numberOfNodesInElement, 0
    );
    this -> integrationPoints.resize(
        globalData.numberOfNodesInElement, vector<double>(2, 0)
    );
    this -> pointsEdges.resize(
        globalData.numberOfNodesInElement * 2, vector<double>(2, 0)
    );
    this -> functionNEdges.resize(
        globalData.numberOfNodesInElement * 2,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> functionsN.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> derivativeDNDKsi.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> derivativeDNDEta.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> H.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> C.resize(
        globalData.numberOfNodesInElement,
        std::vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> derivativeDNDX.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );
    this -> derivativeDNDY.resize(
        globalData.numberOfNodesInElement,
        vector<double>(globalData.numberOfNodesInElement, 0)
    );

    this -> integrationPoints[0][0] = VALUE_MINUS;
    this -> integrationPoints[0][1] = VALUE_MINUS;
    this -> integrationPoints[1][0] = VALUE_PLUS;
    this -> integrationPoints[1][1] = VALUE_MINUS;
    this -> integrationPoints[2][0] = VALUE_PLUS;
    this -> integrationPoints[2][1] = VALUE_PLUS;
    this -> integrationPoints[3][0] = VALUE_MINUS;
    this -> integrationPoints[3][1] = VALUE_PLUS;

    this -> calculateShapeFunctions();
    this -> pointsOnTheEdges();
}


//----------DERIVATIVES OF THE SHAPE FUNCTIONS AFTER KSI AND ETA, SHAPE FUNCTIONS N----------------
void UniversalElement::calculateShapeFunctions() {
    for (int i = 0; i < globalData.numberOfNodesInElement; i++) {
        for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
            shapeFunctionsKsiEta(this -> integrationPoints[i][0], this -> integrationPoints[i][1]);
            derivativeDNDKsi[i][j] = dNdKsi[j];
            derivativeDNDEta[i][j] = dNdEta[j];

            shapeFunctionsN(this -> integrationPoints[i][0], this -> integrationPoints[i][1]);
            functionsN[i][j] = N[j];
        }
    }
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

void UniversalElement::shapeFunctionsN(double ksi, double eta) {
    N[0] = 0.25 * (1 - ksi) * (1 - eta);
    N[1] = 0.25 * (1 + ksi) * (1 - eta);
    N[2] = 0.25 * (1 + ksi) * (1 + eta);
    N[3] = 0.25 * (1 - ksi) * (1 + eta);
}


//-----------------------------POINTS ON THE EDGES-----------------------------
void UniversalElement::pointsOnTheEdges(){

    this -> pointsEdges[0][0] = VALUE_MINUS;
    this -> pointsEdges[0][1] = - 1;
    this -> pointsEdges[1][0] = VALUE_PLUS;
    this -> pointsEdges[1][1] = - 1;
    this -> pointsEdges[2][0] = 1;
    this -> pointsEdges[2][1] = VALUE_MINUS;
    this -> pointsEdges[3][0] = 1;
    this -> pointsEdges[3][1] = VALUE_PLUS;
    this -> pointsEdges[4][0] = VALUE_PLUS;
    this -> pointsEdges[4][1] = 1;
    this -> pointsEdges[5][0] = VALUE_MINUS;
    this -> pointsEdges[5][1] = 1;
    this -> pointsEdges[6][0] = - 1;
    this -> pointsEdges[6][1] = VALUE_PLUS;
    this -> pointsEdges[7][0] = - 1;
    this -> pointsEdges[7][1] = VALUE_MINUS;

    int m = 6;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            shapeFunctionsN(this->pointsEdges[i][0], this->pointsEdges[i][1]);
            functionNEdges[i][j] = N[j];
        }
        m++;
    }
}


//-----------------------------CREATE MATRIX HBC AND VECTOR P-------------------------------------------
void UniversalElement::matrixHBCandVecP(Element elements, Node *node1, Node *node2, double detJ){

    if (elements.nodes[0] -> id == node1 -> id &&
        elements.nodes[globalData.numberOfNodesInElement-1] -> id == node2 -> id) {

        for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
            for (int k = 0; k < globalData.numberOfNodesInElement; k++) {

                HBC[j][k] = ((functionNEdges[globalData.numberOfNodesInElement * 2 - 2][j] *
                              functionNEdges[globalData.numberOfNodesInElement * 2 - 2][k]) +
                             (functionNEdges[globalData.numberOfNodesInElement * 2 - 1][j] *
                              functionNEdges[globalData.numberOfNodesInElement * 2 - 1][k])
                            ) * detJ * globalData.alfa;
            }
            vecP[j] = (((functionNEdges[globalData.numberOfNodesInElement * 2 - 2][j]) +
                         (functionNEdges[globalData.numberOfNodesInElement * 2 - 1][j])) *
                        globalData.ambientTemperature * globalData.alfa * detJ);
        }
    }
    else {
        for (int m = 1; m < globalData.numberOfNodesInElement; m++) {
            if (elements.nodes[m - 1] -> id == node1 -> id &&
                elements.nodes[m] -> id == node2 -> id) {

                for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
                    for (int k = 0; k < globalData.numberOfNodesInElement; k++) {

                        HBC[j][k] += (((functionNEdges[2 * m - 2][j] * functionNEdges[2 * m - 2][k]) +
                                      (functionNEdges[2 * m - 1][j] * functionNEdges[2 * m - 1][k])
                                     ) * detJ * globalData.alfa
                                     );
                    }
                    vecP[j] += (((functionNEdges[2 * m - 2][j]) + (functionNEdges[2 * m - 1][j])) *
                                globalData.ambientTemperature * globalData.alfa * detJ
                                );
                }

            }
        }
    }

}

//--------------------CREATE MATRIXES H AND C----------------------------
void UniversalElement::createMatrixHandC(Element element) {

    for (int i = 0 ; i < 4; i++) {
        for (int j=0; j< 4; j++ ){
            H[i][j] = 0;
            C[i][j] = 0;
            derivativeDNDX[i][j] = 0;
            derivativeDNDY[i][j] = 0;
        }
    }

   for (int i = 0; i < globalData.numberOfNodesInElement; i++) {

       this -> calculcateJacobiTransformation(element, i);

//----------------------DN/DX AND DN/DY-------------------------------
       for (int j = 0; j < globalData.numberOfNodesInElement; j++) {

           derivativeDNDX[i][j] = (
                   jacobiReverseMatrix[0][0] * derivativeDNDKsi[i][j] +
                   jacobiReverseMatrix[0][1] * derivativeDNDEta[i][j]
           );
           derivativeDNDY[i][j] = (
                   jacobiReverseMatrix[1][0] * derivativeDNDKsi[i][j] +
                   jacobiReverseMatrix[1][1] * derivativeDNDEta[i][j]
           );
       }
//-------------------CREATE MATRIXES H AND C-----------------------
       for (int j = 0; j < globalData.numberOfNodesInElement; j++) {
           for (int m = 0; m < globalData.numberOfNodesInElement; m++) {

               H[j][m] += ((derivativeDNDX[i][j] * derivativeDNDX[i][m] +
                            derivativeDNDY[i][j] * derivativeDNDY[i][m]) *
                           globalData.conductivity * det
               );
               C[j][m] += ((functionsN[i][j] * functionsN[i][m]) *
                           globalData.heat * globalData.density * det
               );

           }
       }

   }
}

//-----------------------JACOBI TRANSFORMATION----------------------------
void UniversalElement::calculcateJacobiTransformation(Element element, int p) {
    det = 0;

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++) {
            jacobiMatrix[i][j] = 0;
            jacobiReverseMatrix[i][j] = 0;
        }
    }
//--------------------------CREATE JACOBI MATRIX-------------------------------
    for (int i = 0; i < globalData.numberOfNodesInElement; i++) {

        jacobiMatrix[0][0] += derivativeDNDKsi[p][i] * element.nodes[i]->x;
        jacobiMatrix[0][1] += derivativeDNDKsi[p][i] * element.nodes[i]->y;
        jacobiMatrix[1][0] += derivativeDNDEta[p][i] * element.nodes[i]->x;
        jacobiMatrix[1][1] += derivativeDNDEta[p][i] * element.nodes[i]->y;
    }

//----------------------CREATE DET OF JACOBI MATRIX--------------------------------
    det = jacobiMatrix[0][0] * jacobiMatrix[1][1] - jacobiMatrix[0][1] * jacobiMatrix[1][0];


//----------------------CREATE REVERSE JACOBI MATRIX-------------------------------
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            jacobiReverseMatrix[i][j] = (1 / det) * jacobiMatrix[i][j];
        }
    }
}



//---------------------------------PRINTS-------------------------------------------
void UniversalElement::printLocal() {
    cout << "DERIVATIVE KSI:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << derivativeDNDKsi[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "DERIVATIVES ETA:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << derivativeDNDEta[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "SHAPE FUNCTIONS OF KSI AND ETA:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << functionsN[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "JACOBI MATRIX:" << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cout << jacobiMatrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "DET:" << endl;
    cout << det << endl;

    cout << endl << "REVERSE JACOBI MATRIX:" << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cout << jacobiReverseMatrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "MATRIX H LOCAL:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << H[i][j] << "  ";
        }
        cout << endl;
    }

    cout << endl << "MATRIX C LOCAL:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << C[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    cout << endl << "FUNCTION N VECTOR P:" << endl;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            cout << functionNEdges[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;
}
