//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_GLOBALDATA_H
#define MES_PROJ_GLOBALDATA_H


class GlobalData {
public:
    double width, height;
    int numberOfWidth, numberOfHeight;
    int numberOfElements;
    int numberOfNodes;
    int numberOfNodesInElement;
    double conductivity;
    double heat;
    double density;

    GlobalData();
    void print();


};

#endif //MES_PROJ_GLOBALDATA_H