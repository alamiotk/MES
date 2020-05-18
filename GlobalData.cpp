//
// Created by ala on 08.03.20.

#include "GlobalData.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;


GlobalData::GlobalData(){
    fstream plik;
    plik.open("mes.txt");
    plik >> initialTemperature;
    plik >> simulationTime;
    plik >> simulationStepTime;
    plik >> ambientTemperature;
    plik >> alfa;


    plik >> height;
    plik >> width;
    plik >> numberOfHeight;
    plik >> numberOfWidth;
    plik >> numberOfNodesInElement;
    plik >> conductivity;
    plik >> heat;
    plik >> density;
    this -> numberOfElements = (numberOfHeight-1) * (numberOfWidth-1);
    this -> numberOfNodes = numberOfHeight * numberOfWidth;
}

void GlobalData::print(){
    cout << initialTemperature << endl << simulationTime << endl << simulationStepTime << endl << ambientTemperature << endl;
    cout << alfa << endl << height << endl << width << endl << numberOfHeight << endl << numberOfWidth << endl;
    cout << numberOfNodesInElement << endl << conductivity << endl << heat << endl << density << endl;
}

