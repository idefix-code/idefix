#include <fstream>
#include "idefix.hpp"
#include "input.hpp"
#include "gitversion.h"

Input::Input() {
    std::cout << "Creating Input\n";
    

    /*
    // Number of integrator stages
    nstages=2;

    // First time step
    firstDt = 1.0e-4;

    // Final time
    tfinal = 1.0e-2;

    // Periodicty of outputs
    tperiodVTK=0.1;
    */


}

// Create input from file filename
// This routine expects input file of the following form:
// [Blockname]                                                              # comments2
// Parameter_name     Parametervalue1   Parametervalue2  Parametervalue3... # comments1
//
// Comments are allowed everywhere. Anything after # is ignored in the line
// Blockname should refer to one of Idefix class which will use the parameters in said block
// Everything is stored in a map of maps of vectors of strings :-)

Input::Input(std::string filename) {
    std::ifstream file;
    std::string line, lineWithComments, blockName, paramName, paramValue;
    std::size_t firstChar, lastChar;
    bool haveBlock = false;
    std::stringstream msg;

    try
    {
        file.open(filename);
    }
    catch(...)
    {
        msg << "Input constructor cannot open file " << filename;
        IDEFIX_ERROR(msg);
    }

    while(std::getline(file, lineWithComments)) 
    {
        std::cout << "Interpreting line '"<<lineWithComments <<"'"<<std::endl;
        line = lineWithComments.substr(0, lineWithComments.find("#",0));
        if (line.empty()) continue;         // skip blank line
        firstChar = line.find_first_not_of(" ");
        if (firstChar == std::string::npos) continue;          // line is all white space

        if (line.compare(firstChar, 1, "[") == 0) {              // a new block
            firstChar++;
            lastChar = (line.find_first_of("]", firstChar));

            if (lastChar == std::string::npos) {
                msg << "Block name '" << blockName << "' in file '" << filename << "' not properly ended";
                IDEFIX_ERROR(msg);
            }
            blockName.assign(line, firstChar, lastChar-1);
            haveBlock = true;

            std::cout << "Block " << blockName << std::endl;
            continue;   // Go to next line
        }   // New block

        // At this point, we should have a parameter set in the line
        if(haveBlock == false) {
            msg << "Input file '" << filename << "' must specify a block name before the first parameter";
            IDEFIX_ERROR(msg);
        }
        std::cout << "Parameter:" << line << std::endl;

        std::stringstream streamline(line);
        // Store the name of the parameter
        streamline >> paramName;
        
        // Store the parameters in parameter block
        std::cout << "Pushing input[" << blockName << "][" << paramName << "]:" << std::endl;
        while(streamline >> paramValue) {
            inputParameters[blockName][paramName].push_back(paramValue);
            std::cout << "\t\t [" << inputParameters[blockName][paramName].size() << "]=" << paramValue << std::endl;
        }

    }
    file.close();
}

// This routine prints the parameters stored in the inputParameters structure
void Input::PrintParameters() {
    std::string blockName, paramName, paramValue;
    for(std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.begin(); block != inputParameters.end(); block++ ) {
        blockName=block->first;
        std::cout << "[" << blockName << "]" << std::endl;
        for(std::map<std::string, std::vector<std::string>>::iterator param = block->second.begin(); param !=block->second.end(); param++) {
            paramName=param->first;
            std::cout << "\t" << paramName << "\t";
            for(std::vector<std::string>::iterator value = param->second.begin(); value != param->second.end(); value++){
                paramValue = *value;
                std::cout << "\t" << paramValue;
            } 
            std::cout << std::endl;
        }
    }

}

// Get a string in a block, parameter, position of the file
std::string Input::GetString(std::string blockName, std::string paramName, int num) {
    std::stringstream msg;

    std::string value;

    std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.find(blockName);
    if(block != inputParameters.end()) {
        // Block exists
        std::map<std::string, std::vector<std::string>>::iterator param = block->second.find(paramName);
        if(param != block->second.end()) {
            // Parameter exist
            if(num<param->second.size()) {
                // Vector is long enough
                value=param->second[num];
            }
            else {
                // Vector is not long enough
                msg << "Index " << num << " cannot be found in block:parameter" << blockName << ":" << paramName;
                IDEFIX_ERROR(msg);
            }
        }
        else {
            msg << "Parameter " << paramName << " cannot be found in block [" << blockName <<"]" ;
            IDEFIX_ERROR(msg);
        }
    }
    else {
        msg << "BlockName " << blockName << " cannot be found";
        IDEFIX_ERROR(msg);
    }
    return(value);
}

// Get a real number in a block, parameter, position of the file
real Input::GetReal(std::string blockName, std::string paramName, int num) {
    std::stringstream msg;

    real value;

    std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.find(blockName);
    if(block != inputParameters.end()) {
        // Block exists
        std::map<std::string, std::vector<std::string>>::iterator param = block->second.find(paramName);
        if(param != block->second.end()) {
            // Parameter exist
            if(num<param->second.size()) {
                // Vector is long enough
                #ifdef USE_DOUBLE
                value = std::stod(param->second[num], NULL);
                #else
                value = std::stof(param->second[num], NULL);
                #endif
            }
            else {
                // Vector is not long enough
                msg << "Index " << num << " cannot be found in block:parameter" << blockName << ":" << paramName;
                IDEFIX_ERROR(msg);
            }
        }
        else {
            msg << "Parameter " << paramName << " cannot be found in block [" << blockName <<"]" ;
            IDEFIX_ERROR(msg);
        }
    }
    else {
        msg << "BlockName " << blockName << " cannot be found";
        IDEFIX_ERROR(msg);
    }
    return(value);
}

// Get an integer number in a block, parameter, position of the file
int Input::GetInt(std::string blockName, std::string paramName, int num) {
    std::stringstream msg;

    int value;

    std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.find(blockName);
    if(block != inputParameters.end()) {
        // Block exists
        std::map<std::string, std::vector<std::string>>::iterator param = block->second.find(paramName);
        if(param != block->second.end()) {
            // Parameter exist
            if(num<param->second.size()) {
                // Vector is long enough
                value = std::stoi(param->second[num], NULL);
            }
            else {
                // Vector is not long enough
                msg << "Index " << num << " cannot be found in block:parameter" << blockName << ":" << paramName;
                IDEFIX_ERROR(msg);
            }
        }
        else {
            msg << "Parameter " << paramName << " cannot be found in block [" << blockName <<"]" ;
            IDEFIX_ERROR(msg);
        }
    }
    else {
        msg << "BlockName " << blockName << " cannot be found";
        IDEFIX_ERROR(msg);
    }
    return(value);
}



void Input::PrintLogo() {
    std::cout << "                                  .:HMMMMHn:.  ..:n.."<< std::endl;
    std::cout << "                                .H*'``     `'%HM'''''!x."<< std::endl;
    std::cout << "         :x                    x*`           .(MH:    `#h."<< std::endl;
    std::cout << "        x.`M                   M>        :nMMMMMMMh.     `n."<< std::endl;
    std::cout << "         *kXk..                XL  nnx:.XMMMMMMMMMMML   .. 4X."<< std::endl;
    std::cout << "          )MMMMMx              'M   `^?M*MMMMMMMMMMMM:HMMMHHMM."<< std::endl;
    std::cout << "          MMMMMMMX              ?k    'X ..'*MMMMMMM.#MMMMMMMMMx"<< std::endl;
    std::cout << "         XMMMMMMMX               4:    M:MhHxxHHHx`MMx`MMMMMMMMM>"<< std::endl;
    std::cout << "         XM!`   ?M                `x   4MM'`''``HHhMMX  'MMMMMMMM"<< std::endl;
    std::cout << "         4M      M                 `:   *>     `` .('MX   '*MMMM'"<< std::endl;
    std::cout << "          MX     `X.nnx..                        ..XMx`     'M*X"<< std::endl;
    std::cout << "           ?h.    ''```^'*!Hx.     :Mf     xHMh  M**MMM      4L`"<< std::endl;
    std::cout << "            `*Mx           `'*n.x. 4M>   :M` `` 'M    `       %"<< std::endl;
    std::cout << "             '%                ``*MHMX   X>      !"<< std::endl;
    std::cout << "            :!                    `#MM>  X>      `   :x"<< std::endl;
    std::cout << "           :M                        ?M  `X     .  ..'M"<< std::endl;
    std::cout << "           XX                       .!*X  `x   XM( MMx`h"<< std::endl;
    std::cout << "          'M>::                        `M: `+  MMX XMM `:"<< std::endl;
    std::cout << "          'M> M                         'X    'MMX ?MMk.Xx.."<< std::endl;
    std::cout << "          'M> ?L                     ...:!     MMX.H**'MMMM*h"<< std::endl;
    std::cout << "           M>  #L                  :!'`MM.    . X*`.xHMMMMMnMk."<< std::endl;
    std::cout << "           `!   #h.      :L           XM'*hxHMM*MhHMMMMMMMMMM'#h"<< std::endl;
    std::cout << "           +     XMh:    4!      x   :f   MM'   `*MMMMMMMMMM%  `X"<< std::endl;
    std::cout << "           M     Mf``tHhxHM      M>  4k xxX'      `#MMMMMMMf    `M .>"<< std::endl;
    std::cout << "          :f     M   `MMMMM:     M>   M!MMM:         '*MMf'     'MH*"<< std::endl;
    std::cout << "          !     Xf   'MMMMMX     `X   X>'h.`          :P*Mx.   .d*~.."<< std::endl;
    std::cout << "        :M      X     4MMMMM>     !   X~ `Mh.      .nHL..M#'%nnMhH!'`"<< std::endl;
    std::cout << "       XM      d>     'X`'**h     'h  M   ^'MMHH+*'`  ''''   `'**'"<< std::endl;
    std::cout << "    %nxM>      *x+x.:. XL.. `k     `::X"<< std::endl;
    std::cout << ":nMMHMMM:.  X>  Mn`*MMMMMHM: `:     ?MMn."<< std::endl;
    std::cout << "    `'**MML M>  'MMhMMMMMMMM  #      `M:^*x"<< std::endl;
    std::cout << "         ^*MMttnnMMMMMMMMMMMH>.        M:.4X"<< std::endl;
    std::cout << "                        `MMMM>X   (   .MMM:MM!   ."<< std::endl;
    std::cout << "                          `'''4x.dX  +^ `''MMMMHM?L.."<< std::endl;
    std::cout << "                                ``'           `'`'`'`"<< std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "           This is Idefix " << GITVERSION << std::endl;

}
