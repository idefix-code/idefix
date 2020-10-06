// ********************************************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ********************************************************************************************************

#ifndef INPUT_HPP
#define INPUT_HPP
#include <map>
#include <vector>
#include <iostream>
#include "idefix.hpp"



class Input {
public:
    // Constructor from a file
    Input (int, char ** );
    void PrintParameters();

    // Accessor to input parameters
    std::string GetString(std::string, std::string, int);
    real GetReal(std::string, std::string, int);
    int GetInt(std::string, std::string, int);
    int CheckEntry(std::string, std::string);

    Input();
    void PrintLogo();
private:
    std::string inputFileName;
    std::map<std::string,std::map<std::string,std::vector<std::string>>>   inputParameters;
    void ParseCommandLine(int , char **argv);
};

#endif