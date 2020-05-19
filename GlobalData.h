//
// Created by ala on 06.05.20.
//

#ifndef MES_PROJ_GLOBALDATA_H
#define MES_PROJ_GLOBALDATA_H


class GlobalData {
public:
    double widthGrid, heightGrid;
    unsigned int numberOfNodesWidth, numberOfNodesHeight;
    int numberOfElements;
    int numberOfNodes;
    unsigned int numberOfNodesInElement;
    double conductivity;
    double heat;
    double density;
    int initialTemperature, simulationTime, simulationStepTime;
    int ambientTemperature, alfa;

    GlobalData();
    void print();
};

#endif //MES_PROJ_GLOBALDATA_H
