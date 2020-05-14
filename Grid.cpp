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
            double t = 0;

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




//*****************AGGREGATION H, C AND HBC********************************

void Grid::aggregationHandC() {
    int aggregationMatrixSize =  numberOfHeight * numberOfWidth;

    aggregationMatrixH.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixHBC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));

    UniversalElement universalElement = UniversalElement();

    for (int i = 0; i < numberOfElements; i++) {
      //  cout << endl << "ID: " << i << endl;

        if (this -> checkIfEdge(elements[i])){
            universalElement.matrixHBC(detJ);
        }
        else {
            universalElement.matrixHBC(0);
        }

        universalElement.createMatrixHandC(elements[i]);



        for (int j = 0; j < numberOfNodesInElement; j++){
            for (int k = 0; k < numberOfNodesInElement; k++) {
                int jId = elements[i].nodes[j]->id;
                int kId = elements[i].nodes[k]->id;
                aggregationMatrixH[jId][kId] += universalElement.H[j][k] + universalElement.HBC[j][k];
                aggregationMatrixC[jId][kId] += universalElement.C[j][k];
                aggregationMatrixHBC[jId][kId] += universalElement.HBC[j][k];
            }
        }
    }



    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            aggregationMatrixH[i][j] += aggregationMatrixC[i][j]/50;
           // cout << aggregationMatrixH[i][j] << " ";
        }
      //  cout << endl;
    }
  //  cout << endl;





//********************Print Aggregation H, C AND HBC*****************

    cout << "AGGREGATION: GLOBAL H + HBC + C/DT" << endl;
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
}


bool Grid::checkIfEdge(Element elements){
    if (elements.nodes[0]->borderCondition == 1 &&
        elements.nodes[numberOfNodesInElement-1]->borderCondition == 1) {
        edgeLength(elements.nodes[0], elements.nodes[numberOfNodesInElement - 1]);
    return 1;
    }
    for (int i = 1; i < numberOfNodesInElement; i++) {
        if (elements.nodes[i-1]->borderCondition == 1 &&
            elements.nodes[i]->borderCondition == 1) {
            edgeLength(elements.nodes[i-1], elements.nodes[i]);
            return 1;
        }
    }
    return 0;
}


double Grid::edgeLength(Node *nodes1, Node *nodes2) {
//    double N1X,N1Y,N2X,N2Y;//,N3,N4;
//
//    N1X = 0;
//    N1Y = 0;
//
//    N2X = 0.025;
//    N2Y = 0;

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

