//
// Created by ala on 06.05.20.
//

#include <cmath>
#include "iostream"
#include "Grid.h"
#include "Node.h"
#include "Element.h"
#include "vector"
#include "UniversalElement.h"

using namespace std;

Grid::Grid() {

    this -> deltaX = width / (numberOfWidth - 1);
    this -> deltaY = height / (numberOfHeight - 1);
    this -> detJ = 0;

    int idOfNode = 0;

//*****************CREATES NODES**********************************
    for (int i = 0; i < numberOfWidth; i++) {
        for (int j = 0; j < numberOfHeight; j++, idOfNode++) {

            double x = i * deltaX;
            double y = j * deltaY;
            double t = 100;

            Node newNode = Node(idOfNode, x, y,t,checkBorderCondition(x,y));
            nodes.push_back(newNode);
        }
    }

//****************CREATE ELEMENTS***************************
    int elementsInColumn = numberOfHeight - 1;

    for (int i = 0; i < numberOfElements; i++) {
        int rowElements = i / elementsInColumn * numberOfHeight;
        int id = i % elementsInColumn;

        int a = rowElements + id;
        int b = rowElements + id + numberOfHeight;
        int c = rowElements + id + numberOfHeight + 1;
        int d = rowElements + id + 1;

        vector<Node *> tempNodes;

        tempNodes.push_back(&nodes[a]);
        tempNodes.push_back(&nodes[b]);
        tempNodes.push_back(&nodes[c]);
        tempNodes.push_back(&nodes[d]);

        elements.push_back(Element(tempNodes));

        tempNodes.clear();
    }

}

bool Grid::checkBorderCondition(double x, double y) {
    return (
        x == 0 ||
        y == 0 ||
        x == width ||
        y == height
    );
}




//*****************AGGREGATION H, C, HBC AND VECTOR P, MAIN FUNCTION********************************

void Grid::aggregationHandC() {
    int aggregationMatrixSize =  numberOfHeight * numberOfWidth;

    aggregationMatrixH.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixHBC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationVectoP.resize(aggregationMatrixSize, 0);

    UniversalElement universalElement = UniversalElement();

    for (int i = 0; i < numberOfElements; i++) {
        cout << endl << "ID: " << i << endl;

        for (int g = 0; g < 4; g++){
            for ( int  j = 0; j < 4; j++) {
                universalElement.HBC[g][j] = 0;
                universalElement.vecP[g] = 0;
            }
        }


        double vectorPLocalResult[4] = {};
        double matrixHBLocalResult[4][4] = {};
        this -> checkIfEdge(elements[i], universalElement, vectorPLocalResult, matrixHBLocalResult);

//        cout << "WEKTOR lokalny obliczony: ID:" << i << "; ";
//        for (int f = 0; f < 4; f++) {
//            cout  << vectorPLocalResult[f] << " ";
//        }

        for (int i = 0; i < 4; i++){
            for ( int  j = 0; j < 4; j++) {
                cout << matrixHBLocalResult[i][j] << " ";
            }
            cout << endl;
        }

        universalElement.createMatrixHandC(elements[i]);


        for (int j = 0; j < numberOfNodesInElement; j++){
            for (int k = 0; k < numberOfNodesInElement; k++) {
                int jId = elements[i].nodes[j]->id;
                int kId = elements[i].nodes[k]->id;
                aggregationMatrixH[jId][kId] += universalElement.H[j][k] + matrixHBLocalResult[j][k];
                aggregationMatrixC[jId][kId] += universalElement.C[j][k];
                aggregationMatrixHBC[jId][kId] += matrixHBLocalResult[j][k];
           //     aggregationVectoP[jId] += vectorPLocalResult[j];
                //cout << aggregationVectoP[jId] << endl;
            }
        }
        for ( int i = 0; i < aggregationMatrixSize; i++){
            aggregationVectoP[i] += vectorPLocalResult[i];
        }
    }





    cout << endl << "AGGREGATION: VECTOR P" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        cout << aggregationVectoP[i] << " ";
    }
    cout << endl;


    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            aggregationMatrixH[i][j] += aggregationMatrixC[i][j]/50;
            aggregationVectoP[i] += (aggregationMatrixC[i][j]/50) * 100;
           // cout << aggregationMatrixH[i][j] << " ";
        }
      //  cout << endl;
    }
    cout << endl;


//********************Prints Aggregation H, C AND HBC*****************

    cout << "AGGREGATION: GLOBAL H + HBC + C / DT" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            cout << aggregationMatrixH[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "AGGREGATION: C GLOBAL" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            cout << aggregationMatrixC[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "AGGREGATION: HBC GLOBAL" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            cout << aggregationMatrixHBC[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "AGGREGATION: VECTOR P + C / DT" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        cout << aggregationVectoP[i] << " ";
    }
    cout << endl;



}
void Grid::checkIfEdge(Element elements, UniversalElement universalElement, double vectorPLocalResult[], double matrixHBLocalResult[][4]){
    if (elements.nodes[0]->borderCondition == 1 &&
        elements.nodes[numberOfNodesInElement-1]->borderCondition == 1) {
        edgeLength(elements, elements.nodes[0], elements.nodes[numberOfNodesInElement - 1]);
        universalElement.matrixHBC(elements, elements.nodes[0], elements.nodes[numberOfNodesInElement-1], detJ);
        universalElement.vectorP(elements, elements.nodes[0], elements.nodes[numberOfNodesInElement-1], detJ);
    }

    for (int i = 1; i < numberOfNodesInElement; i++) {
        if (elements.nodes[i-1]->borderCondition == 1 &&
            elements.nodes[i]->borderCondition == 1) {
            edgeLength(elements, elements.nodes[i-1], elements.nodes[i]);
            universalElement.matrixHBC(elements, elements.nodes[i-1], elements.nodes[i], detJ);
            universalElement.vectorP(elements, elements.nodes[i-1], elements.nodes[i], detJ);
        }
    }
    for (int k = 0; k < 4; k++) {
        for (int g = 0; g < 4; g++) {
            vectorPLocalResult[g] = universalElement.vecP[g];
            matrixHBLocalResult[k][g] = universalElement.HBC[k][g];
        }
    }

}


double Grid::edgeLength(Element elements, Node *nodes1, Node *nodes2) {
    double N1X = nodes1 -> x;
    double N1Y = nodes1 -> y;

    double N2X = nodes2 -> x;
    double N2Y = nodes2 -> y;

    detJ = sqrt(pow((N1X-N2X),2)+pow((N1Y-N2Y),2))/2;
    return detJ;
}



//******************PRINTS***********************************

void Grid::printGrid () {

    for(int i = 0; i < numberOfNodes; i ++) {
        cout << " ID: " << nodes[i].id << " X: " << nodes[i].x << " Y: " << nodes[i].y << " t: " << nodes[i].t << " BC: " << nodes[i].borderCondition  << endl << endl;
    }

    cout << "ELEMENTS" << endl;


    for (int i = 0; i < numberOfElements; i++) {
        cout << "ID: " << i << ": " << elements[i].printCoordinates() << endl;
    }
    cout << endl;
}

