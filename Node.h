//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_NODE_H
#define MES_PROJ_NODE_H


class Node {
public:
    int id;
    double x,y,t;
    bool borderCondition;

    Node(int id, double x, double y,double t, bool borderCondition);
};


#endif //MES_PROJ_NODE_H
