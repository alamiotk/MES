//
// Created by ala on 06.05.20.
//

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

    int idOfNode = 0;

    cout << numberOfHeight << " " << numberOfWidth << endl;


//   Create nodes
    for (int i = 0; i < numberOfWidth; i++) {
        for (int j = 0; j < numberOfHeight; j++, idOfNode++) {

            double x = i * deltaX;
            double y = j * deltaY;

            Node newNode = Node(idOfNode, x, y,  checkBorderCondition(x,y));
            nodes.push_back(newNode);

//            cout << newNode.x << " " << newNode.y << " " << newNode.borderCondition << " " << newNode.id << endl;

        }
    }

//    Create elements

    int elementCount = numberOfHeight - 1;
    for (int i = 0; i < numberOfElements; i++) {
        int rowElements = i / elementCount * numberOfHeight;
        int id = i % elementCount;

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
        x == numberOfHeight - 1 ||
        y == numberOfWidth -1
    );
}

void Grid::printGrid () {

    for(int i = 0; i < numberOfNodes; i ++) {
        cout << i << nodes[i].x << " " << nodes[i].y << " " << nodes[i].borderCondition << " " << nodes[i].id << endl;
    }

//    cout << "ELEMENTS" << endl;
//
//
//    for (int i = 0; i < numberOfElements; i++) {
//        cout << i << " " << elements[i].printCoordinates() << endl;
//    }
}

void Grid::aggregationHandC() {
    int aggregationMatrixSize =  numberOfHeight * numberOfWidth;

    aggregationMatrixH.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));
    aggregationMatrixC.resize(aggregationMatrixSize, std::vector<double>(aggregationMatrixSize, 0));


    UniversalElement universalElement = UniversalElement();

//    universalElement.createMatrixHandC(elements[9]);
//    universalElement.print();

    for (int i = 0; i < numberOfElements; i++) {
        universalElement.createMatrixHandC(elements[i]);

        cout << "ID: " << i << endl;
        universalElement.print();

        for (int j = 0; j < numberOfNodesInElement; j++){
            for (int k = 0; k < numberOfNodesInElement; k++) {
                int jId = elements[i].nodes[j]->id;
                int kId = elements[i].nodes[k]->id;
                aggregationMatrixH[jId][kId] += universalElement.H[i][j];
                aggregationMatrixC[jId][kId] += universalElement.C[i][j];
            }
        }
    }


    cout << "macierz H agragcja:" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            cout << aggregationMatrixH[i][j] << " ";
        }
        cout << endl;
    }
//
    cout << "macierz C agragcja:" << endl;
    for (int i = 0; i < aggregationMatrixSize; i++) {
        for (int j = 0; j < aggregationMatrixSize; j++) {
            cout << aggregationMatrixC[i][j] << " ";
        }
        cout << endl;
    }


}