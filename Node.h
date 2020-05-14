//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_NODE_H
#define MES_PROJ_NODE_H


class Node {
public:
//    x,y wynikiem
    int id;
    double x,y,t;
    bool borderCondition;

    Node(int id, double x, double y,double t, bool borderCondition);
//    double getX();
//    double getY();

};


#endif //MES_PROJ_NODE_H
