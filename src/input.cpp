#include <fstream>
#include "idefix.hpp"
#include "input.hpp"
#include "gitversion.h"

Input::Input() {
    idfx::cout << "Creating Input\n";
    

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

Input::Input(std::string filename, int argc, char* argv[] ) {
    std::ifstream file;
    std::string line, lineWithComments, blockName, paramName, paramValue;
    std::size_t firstChar, lastChar;
    bool haveBlock = false;
    std::stringstream msg;

    Input::ParseCommandLine(argc,argv);

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

            continue;   // Go to next line
        }   // New block

        // At this point, we should have a parameter set in the line
        if(haveBlock == false) {
            msg << "Input file '" << filename << "' must specify a block name before the first parameter";
            IDEFIX_ERROR(msg);
        }

        std::stringstream streamline(line);
        // Store the name of the parameter
        streamline >> paramName;
        
        // Store the parameters in parameter block
        while(streamline >> paramValue) {
            inputParameters[blockName][paramName].push_back(paramValue);
        }

    }
    file.close();
}

// This routine parse command line options
void Input::ParseCommandLine(int argc, char **argv) {
    std::stringstream msg;
    for(int i = 1 ; i < argc ; i++) {
        // MPI decomposition argument
        if(std::string(argv[i]) == "-dec") {
            #ifndef WITH_MPI
            IDEFIX_ERROR("Domain decomposition option '-dec' only makes sense when MPI is enabled");
            #endif    
            // Loop on dimensions
            for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
                if ((++i) >= argc) {
                    D_SELECT(msg << "You must specify -dec n1";  ,
                            msg << "You must specify -dec n1 n2";  ,
                            msg << "You must specify -dec n1 n2 n3"; )
                    IDEFIX_ERROR(msg);
                }
                // Store this
                inputParameters["CommandLine"]["dec"].push_back(std::string(argv[i]));
            }
        }
        else {
            msg << "Unknown option " << argv[i];
            IDEFIX_ERROR(msg);
        }
    }
}


// This routine prints the parameters stored in the inputParameters structure
void Input::PrintParameters() {
    std::string blockName, paramName, paramValue;
    idfx::cout << "-----------------------------------------------------------------------------" << std::endl;
    idfx::cout << "Input Parameters:" << std::endl;
    idfx::cout << "-----------------------------------------------------------------------------" << std::endl;
    for(std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.begin(); block != inputParameters.end(); block++ ) {
        blockName=block->first;
        idfx::cout << "[" << blockName << "]" << std::endl;
        for(std::map<std::string, std::vector<std::string>>::iterator param = block->second.begin(); param !=block->second.end(); param++) {
            paramName=param->first;
            idfx::cout << "\t" << paramName << "\t";
            for(std::vector<std::string>::iterator value = param->second.begin(); value != param->second.end(); value++){
                paramValue = *value;
                idfx::cout << "\t" << paramValue;
            } 
            idfx::cout << std::endl;
        }
    }
    idfx::cout << "-----------------------------------------------------------------------------" << std::endl;
    idfx::cout << "-----------------------------------------------------------------------------" << std::endl;

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

// Check that an entry is present in the ini file, if yes, return the number of parameters for given entry
int Input::CheckEntry(std::string blockName, std::string paramName) {
    int result=-1;
    std::map<std::string, std::map<std::string, std::vector<std::string>>>::iterator block = inputParameters.find(blockName);
    if(block != inputParameters.end()) {
        // Block exists
        std::map<std::string, std::vector<std::string>>::iterator param = block->second.find(paramName);
        if(param != block->second.end()) {
            // Parameter exist
            result = param->second.size();
        }
    }
    return(result);
}

void Input::PrintLogo() {
    idfx::cout << "                                  .:HMMMMHn:.  ..:n.."<< std::endl;
    idfx::cout << "                                .H*'``     `'%HM'''''!x."<< std::endl;
    idfx::cout << "         :x                    x*`           .(MH:    `#h."<< std::endl;
    idfx::cout << "        x.`M                   M>        :nMMMMMMMh.     `n."<< std::endl;
    idfx::cout << "         *kXk..                XL  nnx:.XMMMMMMMMMMML   .. 4X."<< std::endl;
    idfx::cout << "          )MMMMMx              'M   `^?M*MMMMMMMMMMMM:HMMMHHMM."<< std::endl;
    idfx::cout << "          MMMMMMMX              ?k    'X ..'*MMMMMMM.#MMMMMMMMMx"<< std::endl;
    idfx::cout << "         XMMMMMMMX               4:    M:MhHxxHHHx`MMx`MMMMMMMMM>"<< std::endl;
    idfx::cout << "         XM!`   ?M                `x   4MM'`''``HHhMMX  'MMMMMMMM"<< std::endl;
    idfx::cout << "         4M      M                 `:   *>     `` .('MX   '*MMMM'"<< std::endl;
    idfx::cout << "          MX     `X.nnx..                        ..XMx`     'M*X"<< std::endl;
    idfx::cout << "           ?h.    ''```^'*!Hx.     :Mf     xHMh  M**MMM      4L`"<< std::endl;
    idfx::cout << "            `*Mx           `'*n.x. 4M>   :M` `` 'M    `       %"<< std::endl;
    idfx::cout << "             '%                ``*MHMX   X>      !"<< std::endl;
    idfx::cout << "            :!                    `#MM>  X>      `   :x"<< std::endl;
    idfx::cout << "           :M                        ?M  `X     .  ..'M"<< std::endl;
    idfx::cout << "           XX                       .!*X  `x   XM( MMx`h"<< std::endl;
    idfx::cout << "          'M>::                        `M: `+  MMX XMM `:"<< std::endl;
    idfx::cout << "          'M> M                         'X    'MMX ?MMk.Xx.."<< std::endl;
    idfx::cout << "          'M> ?L                     ...:!     MMX.H**'MMMM*h"<< std::endl;
    idfx::cout << "           M>  #L                  :!'`MM.    . X*`.xHMMMMMnMk."<< std::endl;
    idfx::cout << "           `!   #h.      :L           XM'*hxHMM*MhHMMMMMMMMMM'#h"<< std::endl;
    idfx::cout << "           +     XMh:    4!      x   :f   MM'   `*MMMMMMMMMM%  `X"<< std::endl;
    idfx::cout << "           M     Mf``tHhxHM      M>  4k xxX'      `#MMMMMMMf    `M .>"<< std::endl;
    idfx::cout << "          :f     M   `MMMMM:     M>   M!MMM:         '*MMf'     'MH*"<< std::endl;
    idfx::cout << "          !     Xf   'MMMMMX     `X   X>'h.`          :P*Mx.   .d*~.."<< std::endl;
    idfx::cout << "        :M      X     4MMMMM>     !   X~ `Mh.      .nHL..M#'%nnMhH!'`"<< std::endl;
    idfx::cout << "       XM      d>     'X`'**h     'h  M   ^'MMHH+*'`  ''''   `'**'"<< std::endl;
    idfx::cout << "    %nxM>      *x+x.:. XL.. `k     `::X"<< std::endl;
    idfx::cout << ":nMMHMMM:.  X>  Mn`*MMMMMHM: `:     ?MMn."<< std::endl;
    idfx::cout << "    `'**MML M>  'MMhMMMMMMMM  #      `M:^*x"<< std::endl;
    idfx::cout << "         ^*MMttnnMMMMMMMMMMMH>.        M:.4X"<< std::endl;
    idfx::cout << "                        `MMMM>X   (   .MMM:MM!   ."<< std::endl;
    idfx::cout << "                          `'''4x.dX  +^ `''MMMMHM?L.."<< std::endl;
    idfx::cout << "                                ``'           `'`'`'`"<< std::endl;
    idfx::cout << std::endl;
    idfx::cout << std::endl;
    idfx::cout << std::endl;
    idfx::cout << "           This is Idefix " << GITVERSION << std::endl;
#ifdef KOKKOS_ENABLE_CUDA
    idfx::cout << "               Compiled for GPU (nvidia-CUDA) " << std::endl;
#else
    idfx::cout << "               Compiled for CPUs " << std::endl;
#endif

}
