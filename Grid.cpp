//
// Created by ala on 06.05.20.
//

#include <cmath>
#include <iomanip>
#include "iostream"
#include "Grid.h"
#include "Node.h"
#include "Element.h"
#include "vector"
#include "UniversalElement.h"

using namespace std;

const double EPS = 1e-10;

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

//********************CHECKING BORDER CONDITIONS TO CREATE NODES****************************
bool Grid::checkBorderCondition(double x, double y) {
    return (
        x == 0 ||
        y == 0 ||
        x == width ||
        y == height
    );
}


//*****************AGGREGATION H, C, HBC AND VECTOR P, MAIN FUNCTION********************************
void Grid::aggregation() {
    int aggregationMatrixSize =  numberOfHeight * numberOfWidth;

    aggregationMatrixH.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixHBC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    matrixCxT0.resize(aggregationMatrixSize, 0);
    aggregationVectoP.resize(aggregationMatrixSize, 0);
    temperatureInitialMatrix.resize(aggregationMatrixSize, initialTemperature);
    temperatureT1.resize(aggregationMatrixSize, 0);
    double minimalTemperature = 0;
    double maximalTemperature = 0;

    UniversalElement universalElement = UniversalElement();



    for (int i = 0; i < numberOfElements; i++) {

        // FILL MATRIX HBC AND VECP WITH ZEROS
        for (int g = 0; g < numberOfNodesInElement; g++){
            for ( int  j = 0; j < numberOfNodesInElement; j++) {
                universalElement.HBC[g][j] = 0;
            }
            universalElement.vecP[g] = 0;
        }

        // CREATE LOCAL VEC AND MATRIX FOR THE RESULTS OF BORDER CONDITIONS
        double vectorPLocalResult[4] = {};
        double matrixHBLocalResult[4][4] = {};
        this -> checkIfEdge(elements[i], universalElement, vectorPLocalResult, matrixHBLocalResult);


        // CREATE MATRIXES H AND C
        universalElement.createMatrixHandC(elements[i]);


        // AGGREGATION OF H, HBC, C, VEC P
        for (int j = 0; j < numberOfNodesInElement; j++){
            int jId = elements[i].nodes[j]->id;
            for (int k = 0; k < numberOfNodesInElement; k++) {

                int kId = elements[i].nodes[k]->id;
                aggregationMatrixH[jId][kId] += universalElement.H[j][k] + matrixHBLocalResult[j][k];
                aggregationMatrixC[jId][kId] += universalElement.C[j][k];
                aggregationMatrixHBC[jId][kId] += matrixHBLocalResult[j][k];
            }
            aggregationVectoP[jId] += vectorPLocalResult[j];
        }

    }

    //************************H (with HBC) + (C/DT)*********************************
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            aggregationMatrixH[i][j] += aggregationMatrixC[i][j]/simulationStepTime;
        }
    }


    //*************************MATRIX (C/DT) x T0***************************************
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            matrixCxT0[i] += (aggregationMatrixC[i][j]/simulationStepTime) * temperatureInitialMatrix[j];
        }
    }

    //***********************VEC P + (C/DT) x T0*****************************************
    for ( int i = 0; i < aggregationMatrixSize; i++) {
        aggregationVectoP[i] += matrixCxT0[i];
    }


    //***********************SOLVING THE EQUATION**********************************************
    this -> solveEquation(aggregationMatrixH, aggregationVectoP, aggregationMatrixSize, temperatureT1);

    for (int b = 0; b < aggregationMatrixSize; b++){
        cout << temperatureT1[b] << " ";
    }
    cout << endl;
    this -> min(temperatureT1, aggregationMatrixSize, minimalTemperature);
    this -> max(temperatureT1, aggregationMatrixSize, maximalTemperature);

    cout << maximalTemperature << " " << minimalTemperature << endl;

//********************Prints Aggregation H, C AND HBC*****************

//    cout << "AGGREGATION: GLOBAL H + HBC + C / DT" << endl;
//    for (int i = 0; i < aggregationMatrixSize; i++) {
//        for (int j = 0; j < aggregationMatrixSize; j++) {
//            cout << aggregationMatrixH[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
//
//    cout << "AGGREGATION: C GLOBAL" << endl;
//    for (int i = 0; i < aggregationMatrixSize; i++) {
//        for (int j = 0; j < aggregationMatrixSize; j++) {
//            cout << aggregationMatrixC[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    cout << "AGGREGATION: HBC GLOBAL" << endl;
//    for (int i = 0; i < aggregationMatrixSize; i++) {
//        for (int j = 0; j < aggregationMatrixSize; j++) {
//            cout << aggregationMatrixHBC[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    cout << "AGGREGATION: VECTOR P + (C / DT)*T0" << endl;
//    for (int i = 0; i < aggregationMatrixSize; i++) {
//        cout << aggregationVectoP[i] << " ";
//    }
//    cout << endl;

}
double Grid::min(vector<double> temperatureT1, int aggregationMatrixSize, double minimalTemperature) {
    double min = 10000;
    for (int i = 0; i < aggregationMatrixSize; i++){
        if (temperatureT1[i] < min) {
            min = temperatureT1[i];
        }
    }
    minimalTemperature = min;
    return minimalTemperature;

}
double Grid::max(vector<double> temperatureT1, int aggregationMatrixSize,  double maximalTemperature) {
    double max = 0;
    for (int i = 0; i < aggregationMatrixSize; i++){
        if (temperatureT1[i] > max) {
            max = temperatureT1[i];
        }
    }
    maximalTemperature = max;
    return maximalTemperature;
}



//*************SOLUTION OF THE EQUATION SYSTEM BY THE GAUSS ELIMINATION METHOD***********************
vector<double> Grid::solveEquation(vector<vector<double> > aggregationMatrixH, vector<double> aggregationVectoP, int aggregationMatrixSize, vector<double> temperatureT1) {
    for (int p = 0; p < aggregationMatrixSize; p++){
        // FIND PIVOT ROW AND SWAP
        int max = p;
        for ( int i = p + 1; i < aggregationMatrixSize; i++) {
            if (fabs(aggregationMatrixH[i][p]) > fabs(aggregationMatrixH[max][p])) {
                max = i;
            }
        }

        // PIVOT WITHIN AGGREGATION H AND AGGREGATIONVECTOR P -> A & B
        vector<double> temp = aggregationMatrixH[p];
        aggregationMatrixH[p] = aggregationMatrixH[max];
        aggregationMatrixH[max] = temp;

        double t = aggregationVectoP[p];
        aggregationVectoP[p] = aggregationVectoP[max];
        aggregationVectoP[max] = t;

        if ( fabs(aggregationMatrixH[p][p]) <= EPS) {
            break;
        }

        for ( int i = p + 1; i < aggregationMatrixSize; i++) {
            double alpha = aggregationMatrixH[i][p] / aggregationMatrixH[p][p];
            aggregationVectoP[i] -= alpha * aggregationVectoP[p];
            for ( int j = p; j < aggregationMatrixSize; j++){
                aggregationMatrixH[i][j] -= alpha * aggregationMatrixH[p][j];
            }
        }
    }

    // BACK SUBSTITUTION
    vector<double> x;
    x.resize(aggregationMatrixSize,0);

    for(int i = aggregationMatrixSize - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < aggregationMatrixSize; j++){
            sum += aggregationMatrixH[i][j] * x[j];
        }
        x[i] = (aggregationVectoP[i] - sum) / aggregationMatrixH[i][i];
    }

    for (int i = 0; i < aggregationMatrixSize; i++){
       temperatureT1[i] = x[i];
    }
    return temperatureT1;
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
    for (int k = 0; k < numberOfNodesInElement; k++) {
        for (int g = 0; g < numberOfNodesInElement; g++) {
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

