//
// Created by ala on 06.05.20.
//

#include "Node.h"

Node::Node(int id, double x, double y,double t, bool borderCondition) {
    this -> x = x;
    this -> y = y;
    this -> id = id;
    this -> t = t;
    this -> borderCondition = borderCondition;
}
