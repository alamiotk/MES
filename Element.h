//
// Created by ala on 06.05.20.
//

#include <vector>
#include "Node.h"
#include <string>


#ifndef MES_PROJ_ELEMENT_H
#define MES_PROJ_ELEMENT_H

using namespace std;

class Element {
public:
    vector<Node *> nodes;

    Element(vector<Node *> nodes);
    string printNodesInElement();
};


#endif //MES_PROJ_ELEMENT_H
