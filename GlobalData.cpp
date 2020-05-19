//
// Created by ala on 08.03.20.

#include "GlobalData.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;


GlobalData::GlobalData(){

//-------------LOADING DATA-------------------------
    fstream plik;
    plik.open("mes1.txt");
    //plik.open("mes2.txt");
    plik >> initialTemperature;
    plik >> simulationTime;
    plik >> simulationStepTime;
    plik >> ambientTemperature;
    plik >> alfa;
    plik >> heightGrid;
    plik >> widthGrid;
    plik >> numberOfNodesHeight;
    plik >> numberOfNodesWidth;
    plik >> numberOfNodesInElement;
    plik >> conductivity;
    plik >> heat;
    plik >> density;

    this -> numberOfElements = (numberOfNodesHeight-1) * (numberOfNodesWidth-1);
    this -> numberOfNodes = numberOfNodesHeight * numberOfNodesWidth;
}

//-------------PRINT DATA--------------------------
void GlobalData::print(){
    cout << initialTemperature << endl << simulationTime << endl << simulationStepTime << endl << ambientTemperature << endl;
    cout << alfa << endl << heightGrid << endl << widthGrid << endl << numberOfNodesHeight << endl << numberOfNodesWidth << endl;
    cout << numberOfNodesInElement << endl << conductivity << endl << heat << endl << density << endl;
}

