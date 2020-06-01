//
// Created by ala on 06.05.20.
//

#include "iostream"
#include "vector"

#include "Node.h"
#include "Element.h"

using namespace std;

Element::Element(vector<Node *> nodes) {
    this -> nodes = nodes;
}

//-----------PRINT NODES IN ELEMENT------------------

string Element::printNodesInElement() {
    string s;
    for (int i = 0; i < nodes.size(); i++) {
        s += to_string(nodes[i] -> id) + " ";
    }
    return s;
}