#ifndef INPUT_HPP
#define INPUT_HPP
#include <map>
#include <vector>
#include <iostream>
#include "idefix.hpp"



class Input {
public:
    // Constructor from a file
    Input (std::string );
    void PrintParameters();

    // Accessor to input parameters
    std::string GetString(std::string, std::string, int);
    real GetReal(std::string, std::string, int);
    int GetInt(std::string, std::string, int);
    /* 
    
    int npoints[3];
    real xstart[3];
    real xend[3];
    int nghost[3];
    int nstages;
    real firstDt;
    real tperiodVTK;
    real tfinal;
    */
    Input();
    void PrintLogo();
private:
    
    std::map<std::string,std::map<std::string,std::vector<std::string>>>   inputParameters;
};

#endif