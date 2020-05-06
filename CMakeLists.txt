cmake_minimum_required(VERSION 3.12)
project(Mes_Proj)

set(CMAKE_CXX_STANDARD 14)

add_executable(Mes_Proj main.cpp GlobalData.cpp GlobalData.h Node.cpp Node.h Grid.cpp Grid.h Element.cpp Element.h UniversalElement.cpp UniversalElement.h)