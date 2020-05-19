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
//------------SET WIDTH AND HEIGHT OF EACH ELEMENT-----------------------
    this -> deltaX = widthGrid / (numberOfNodesWidth - 1);
    this -> deltaY = heightGrid / (numberOfNodesHeight - 1);
    this -> detJ = 0;

//----------------------CREATE NODES--------------------------------------
    int idOfNode = 0;
    for (int i = 0; i < numberOfNodesWidth; i++) {
        for (int j = 0; j < numberOfNodesHeight; j++, idOfNode++) {
            double x = i * deltaX;
            double y = j * deltaY;
            double t = 0;

            Node newNode = Node(idOfNode, x, y,t,checkBorderCondition(x,y));
            nodes.push_back(newNode);
        }
    }

//---------------------CREATE ELEMENTS----------------------------------------
    int elementsInColumn = numberOfNodesHeight - 1;

    for (int i = 0; i < numberOfElements; i++) {
        int rowElements = i / elementsInColumn * numberOfNodesHeight;
        int id = i % elementsInColumn;

        int a = rowElements + id;
        int b = rowElements + id + numberOfNodesHeight;
        int c = rowElements + id + numberOfNodesHeight + 1;
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

//----------------CHECKING BORDER CONDITIONS TO CREATE NODES----------------------
bool Grid::checkBorderCondition(double x, double y) {
    return (
        x == 0 ||
        y == 0 ||
        x == widthGrid ||
        y == heightGrid
    );
}


//----------------AGGREGATION H, C, HBC AND VECTOR P, SOLVE EQUATION-----------------
void Grid::aggregation() {
    aggregationMatrixSize =  numberOfNodesHeight * numberOfNodesWidth;
    aggregationMatrixH.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationVectoP.resize(aggregationMatrixSize, 0);

    temperatureInitialMatrix.resize(aggregationMatrixSize, initialTemperature);
    matrixCxT0.resize(aggregationMatrixSize, 0);
    temperatureT1.resize(aggregationMatrixSize, 0);


    double minimalTemperature;
    double maximalTemperature;

    UniversalElement universalElement = UniversalElement();

//------------------LOOP BY TIME----------------------------
    for(int r = 0; r < (simulationTime / simulationStepTime); r++) {

//---------FILL ALL AGGREGATION MATRIXES WITH ZEROS---------
        for (int q = 0; q < aggregationMatrixSize; q++) {
            for (int y = 0; y < aggregationMatrixSize; y++) {
               aggregationMatrixH[q][y] = 0;
               aggregationMatrixC[q][y] = 0;
            }
            aggregationVectoP[q] = 0;
            matrixCxT0[q] = 0;
        }


//----------------LOOP BY ELEMENTS-----------------
        for (int i = 0; i < numberOfElements; i++) {

//------------CREATE LOCAL VEC AND MATRIX FOR THE RESULTS OF BORDER CONDITIONS----------
            double vectorPLocalResult[4] = {};
            double matrixHBLocalResult[4][4] = {};

            this -> checkIfEdge(elements[i], universalElement, vectorPLocalResult, matrixHBLocalResult);


//--------------- CREATE LOCAL MATRIXES H AND C--------------------------
            universalElement.createMatrixHandC(elements[i]);


//----------------AGGREGATION OF H, HBC, C, VEC P-----------------------------------
            for (int j = 0; j < numberOfNodesInElement; j++) {
                int jId = elements[i].nodes[j]->id;
                for (int k = 0; k < numberOfNodesInElement; k++) {
                    int kId = elements[i].nodes[k]->id;
                    aggregationMatrixH[jId][kId] += universalElement.H[j][k] + matrixHBLocalResult[j][k];
                    aggregationMatrixC[jId][kId] += universalElement.C[j][k];
                }
                aggregationVectoP[jId] += vectorPLocalResult[j];
            }
        }

//-----------------------H (with HBC) + (C/DT)--------------------------
//------------------------MATRIX (C/DT) x T0----------------------------
//------------------------VEC P + (C/DT) x T0---------------------------
        for (int i = 0; i < aggregationMatrixSize; i++) {
            for (int j = 0; j < aggregationMatrixSize; j++) {
                aggregationMatrixH[i][j] += aggregationMatrixC[i][j] / simulationStepTime;
                matrixCxT0[i] += (aggregationMatrixC[i][j] / simulationStepTime) *
                temperatureInitialMatrix[j];
            }
            aggregationVectoP[i] += matrixCxT0[i];
        }

        //this -> printVectorP();


//----------------------SOLVING THE EQUATION-------------------------------------
        temperatureT1 = this -> solveEquation(aggregationMatrixH, aggregationVectoP, aggregationMatrixSize);

        minimalTemperature = this -> min(temperatureT1, aggregationMatrixSize);
        maximalTemperature = this -> max(temperatureT1, aggregationMatrixSize);

        cout << "ITERACJA: " << r << endl;
        this -> printTemperatureT1(minimalTemperature, maximalTemperature);


        for (int i = 0; i < aggregationMatrixSize; i++){
            nodes[i].t = temperatureT1[i];
        }

        //this -> printGrid();

//------------CHANGE T0 TO T1, FILL T1 WITH ZEROS-------------------------
        for (int i = 0; i < aggregationMatrixSize; i++) {
            temperatureInitialMatrix[i] = temperatureT1[i];
            temperatureT1[i] = 0;
        }
        //universalElement.printLocal();
    }
    //this -> aggregationPrint();
}


//---------------------CHECK BORDER CONDITIONS AND CALCULATE VEC P AND HBC----------------------------
void Grid::checkIfEdge(Element elements, UniversalElement universalElement, double vectorPLocalResult[], double matrixHBLocalResult[][4]){
    if (elements.nodes[0]->borderCondition == 1 &&
        elements.nodes[numberOfNodesInElement-1]->borderCondition == 1) {

        edgeLength(elements, elements.nodes[0], elements.nodes[numberOfNodesInElement - 1]);
        universalElement.matrixHBCandVecP(elements, elements.nodes[0], elements.nodes[numberOfNodesInElement - 1], detJ);
    }

    for (int i = 1; i < numberOfNodesInElement; i++) {
        if (elements.nodes[i-1]->borderCondition == 1 &&
            elements.nodes[i]->borderCondition == 1) {

            edgeLength(elements, elements.nodes[i-1], elements.nodes[i]);
            universalElement.matrixHBCandVecP(elements, elements.nodes[i - 1], elements.nodes[i], detJ);
        }
    }

//---------------------ASSIGNMENT-------------------
    for (int k = 0; k < numberOfNodesInElement; k++) {
        for (int g = 0; g < numberOfNodesInElement; g++) {
            vectorPLocalResult[g] = universalElement.vecP[g];
            matrixHBLocalResult[k][g] = universalElement.HBC[k][g];
        }
    }
}

//------------------------CALCULATE DET J-------------------------------
double Grid::edgeLength(Element elements, Node *nodes1, Node *nodes2) {
    double N1X = nodes1 -> x;
    double N1Y = nodes1 -> y;

    double N2X = nodes2 -> x;
    double N2Y = nodes2 -> y;

    detJ = sqrt(pow((N1X - N2X), 2) + pow((N1Y - N2Y), 2)) / 2;
    return detJ;
}




//----------------SOLVING THE EQUATION BY THE GAUSS ELIMINATION METHOD---------------------
vector<double> Grid::solveEquation(vector<vector<double> > aggregationMatrixH, vector<double> aggregationVectoP, int aggregationMatrixSize) {

    for (int p = 0; p < aggregationMatrixSize; p++){

//--------------FIND PIVOT ROW AND SWAP-------------------
        int max = p;
        for ( int i = p + 1; i < aggregationMatrixSize; i++) {
            if (fabs(aggregationMatrixH[i][p]) > fabs(aggregationMatrixH[max][p])) {
                max = i;
            }
        }

//------------PIVOT WITHIN AGGREGATION H AND AGGREGATION VECTOR P -> A & B---------
        vector<double> temp = aggregationMatrixH[p];
        aggregationMatrixH[p] = aggregationMatrixH[max];
        aggregationMatrixH[max] = temp;

        double t = aggregationVectoP[p];
        aggregationVectoP[p] = aggregationVectoP[max];
        aggregationVectoP[max] = t;

        if ( fabs(aggregationMatrixH[p][p]) <= EPS) {
            cout << "error" << endl;
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

//-------------BACK SUBSTITUTION------------------------
    vector<double> x;
    x.resize(aggregationMatrixSize,0);

    for(int i = aggregationMatrixSize - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < aggregationMatrixSize; j++){
            sum += aggregationMatrixH[i][j] * x[j];
        }
        x[i] = (aggregationVectoP[i] - sum) / aggregationMatrixH[i][i];
    }

    return x;
}


//---------------------MIN AND MAX TEMPERATURE-----------------------
double Grid::min(vector<double> temperatureT1, int aggregationMatrixSize) {
    double min = 1000000;
    for (int i = 0; i < aggregationMatrixSize; i++){
        if (temperatureT1[i] < min) {
            min = temperatureT1[i];
        }
    }
    return min;

}
double Grid::max(vector<double> temperatureT1, int aggregationMatrixSize) {
    double max = 0;
    for (int i = 0; i < aggregationMatrixSize; i++){
        if (temperatureT1[i] > max) {
            max = temperatureT1[i];
        }
    }
    return max;
}



//------------------------PRINTS-------------------------------------
void Grid::printGrid () {
    cout << "ALL NODES WITH ATRIBUTES:" << endl;
    for(int i = 0; i < numberOfNodes; i ++) {
        cout << "ID: " << nodes[i].id << " X: " << nodes[i].x << " Y: " << nodes[i].y << " t: " << nodes[i].t << " BC: " << nodes[i].borderCondition  << endl;
    }
    cout << endl;

    cout << "ELEMENTS:" << endl;
    for (int i = 0; i < numberOfElements; i++) {
        cout << "ID: " << i << ": " << elements[i].printNodesInElement() << endl;
    }
    cout << endl;
}


void Grid::aggregationPrint() {
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
}

void Grid::printVectorP() {
    cout << endl << "AGGREGATION: VECTOR P + (C / DT)*T0" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        cout << aggregationVectoP[i] << " ";
    }
    cout << endl << endl;
}

void Grid::printTemperatureT1(double minimalTemperature, double maximalTemperature) {
//    cout << "TEMPERATURE T1: ";
//    for (int b = 0; b < aggregationMatrixSize; b++) {
//        cout << temperatureT1[b] << " ";
//    }
//    cout << endl;

    cout << "MIN: " << minimalTemperature << endl << "MAX: " << maximalTemperature << endl << endl;
}