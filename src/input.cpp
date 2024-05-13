// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <dirent.h>

#include <fstream>
#include <string>
#include <csignal>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>
#include <memory>

#include "idefix.hpp"
#include "input.hpp"
#include "version.hpp"
#include "profiler.hpp"

// Flag will be set if a signal has been received
bool Input::abortRequested = false;

Input::Input() {
}

// Create input from file filename
// This routine expects input file of the following form:
// [Blockname]                                # comments2
// Parameter_name   Parametervalue1   Parametervalue2  Parametervalue3... # comments1
//
// Comments are allowed everywhere. Anything after # is ignored in the line
// Blockname should refer to one of Idefix class which will use the parameters
// in said block. Everything is stored in a map of maps of vectors of strings :-)

Input::Input(int argc, char* argv[] ) {
  std::ifstream file;
  std::string line, lineWithComments, blockName, paramName, paramValue;
  std::size_t firstChar, lastChar;
  bool haveBlock = false;
  std::stringstream msg;
  int nParameters = 0;    // # of parameters in current block

  // Tell the system we want to catch the SIGUSR2 signals
  signal(SIGUSR2, signalHandler);

  // Tell the time when input was initialised
  timer.reset();
  lastStopFileCheck = timer.seconds();

  // Default input file name
  this->inputFileName = std::string("idefix.ini");

  // Parse command line (may replace the input file)
  ParseCommandLine(argc,argv);

  file.open(this->inputFileName);

  if(!file) {
    msg << "Input: cannot open input file " << this->inputFileName;
    IDEFIX_ERROR(msg);
  }

  while(std::getline(file, lineWithComments)) {
    line = lineWithComments.substr(0, lineWithComments.find("#",0));
    if (line.empty()) continue;     // skip blank line
    firstChar = line.find_first_not_of(" ");
    if (firstChar == std::string::npos) continue;      // line is all white space

    if (line.compare(firstChar, 1, "[") == 0) {        // a new block
      firstChar++;
      lastChar = (line.find_first_of("]", firstChar));

      if (lastChar == std::string::npos) {
        msg << "Block name '" << blockName << "' in file '"
            << this->inputFileName << "' not properly ended";
        IDEFIX_ERROR(msg);
      }
      // Check if previous block was empty
      if(haveBlock && nParameters == 0) {
        IDEFIX_WARNING(blockName+std::string(" block is empty in "+inputFileName));
        inputParameters[blockName]["!!empty!!"].push_back("!!empty!!");
      }
      blockName.assign(line, firstChar, lastChar-1);
      haveBlock = true;
      nParameters = 0;

      continue;   // Go to next line
    }   // New block

    // At this point, we should have a parameter set in the line
    if(haveBlock == false) {
      msg << "Input file '" << this->inputFileName
          << "' must specify a block name before the first parameter";
      IDEFIX_ERROR(msg);
    }

    std::stringstream streamline(line);
    // Store the name of the parameter
    streamline >> paramName;
    nParameters++;
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
  bool enableLogs = true;
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
    } else if(std::string(argv[i]) == "-restart") {
      std::string sirestart{};
      bool explicitDump = true;     // by default, assume a restart dump # was given
      // Check whether -restart was given with a number or not

      // -restart was the very last parameter
      if((i+1)>= argc) {
        explicitDump = false;
      } else if(std::isdigit(argv[i+1][0]) == 0) {
        // next argiment is another parameter (does not start with a number)
        explicitDump = false;
      }

      if(explicitDump) {
        sirestart = std::string(argv[++i]);
      } else {
        // implicitly restart from the latest existing dumpfile
        sirestart = "-1";
      }
      inputParameters["CommandLine"]["restart"].push_back(sirestart);
      this->restartRequested = true;
      this->restartFileNumber = std::stoi(sirestart);
    } else if(std::string(argv[i]) == "-i") {
      // Loop on dimensions
      if((++i) >= argc) IDEFIX_ERROR(
                      "You must specify -i filename where filename is the name of the input file.");
      this->inputFileName = std::string(argv[i]);
    } else if(std::string(argv[i]) == "-maxcycles") {
      if((i+1)>= argc) {
        IDEFIX_ERROR("-maxcycles requires an additional integer parameter");
      } else if(std::isdigit(argv[i+1][0]) == 0) {
        // next argiment is another parameter (does not start with a number)
        IDEFIX_ERROR("-maxcycles requires an additional integer paramater");
      }
      this->maxCycles = std::stoi(std::string(argv[++i]));
      inputParameters["CommandLine"]["maxCycles"].push_back(std::to_string(maxCycles));
    } else if(std::string(argv[i]) == "-force_init") {
      this->forceInitRequested = true;
    } else if(std::string(argv[i]) == "-nowrite") {
      this->forceNoWrite = true;
      enableLogs = false;
    } else if(std::string(argv[i]) == "-nolog") {
      enableLogs = false;
    } else if(std::string(argv[i]) == "-profile") {
      idfx::prof.EnablePerformanceProfiling();
    } else if(std::string(argv[i]) == "-Werror") {
      idfx::warningsAreErrors = true;
    } else if(std::string(argv[i]) == "-version" || std::string(argv[i]) == "-v") {
      PrintVersion();
      idfx::safeExit(0);
    } else if(std::string(argv[i]) == "-help" || std::string(argv[i]) == "-h") {
      PrintOptions();
      idfx::safeExit(0);
    } else {
      PrintOptions();
      msg << "Unknown option " << argv[i];
      IDEFIX_ERROR(msg);
    }
  }
  if(enableLogs) {
    idfx::cout.enableLogFile();
  }
}


// This routine prints the parameters stored in the inputParameters structure
void Input::ShowConfig() {
  std::string blockName, paramName, paramValue;
  idfx::cout << "-----------------------------------------------------------------------------"
             << std::endl;
  idfx::cout << "Input Parameters using input file " << this->inputFileName << ":" << std::endl;
  idfx::cout << "-----------------------------------------------------------------------------"
             << std::endl;
  for(IdefixInputContainer::iterator block = inputParameters.begin();
      block != inputParameters.end();
      block++ ) {
    blockName=block->first;
    idfx::cout << "[" << blockName << "]" << std::endl;
    for(IdefixBlockContainer::iterator param = block->second.begin();
          param !=block->second.end(); param++) {
      paramName=param->first;
      idfx::cout << "\t" << paramName << "\t";
      for(IdefixParamContainer::iterator value = param->second.begin();
          value != param->second.end(); value++) {
        paramValue = *value;
        idfx::cout << "\t" << paramValue;
      }
      idfx::cout << std::endl;
    }
  }
  idfx::cout << "-----------------------------------------------------------------------------"
             << std::endl;
  idfx::cout << "-----------------------------------------------------------------------------"
             << std::endl;

  #ifdef SINGLE_PRECISION
    idfx::cout << "Input: Compiled with SINGLE PRECISION arithmetic." << std::endl;
  #else
    idfx::cout << "Input: Compiled with DOUBLE PRECISION arithmetic." << std::endl;
  #endif
  // Show dimensionality and other general options:
  idfx::cout << "Input: DIMENSIONS=" << DIMENSIONS << "." << std::endl;
  idfx::cout << "Input: COMPONENTS=" << COMPONENTS << "." << std::endl;
  #ifdef WITH_MPI
    idfx::cout << "Input: MPI ENABLED." << std::endl;
  #endif
  #ifdef KOKKOS_ENABLE_HIP
    idfx::cout << "Input: Kokkos HIP target ENABLED." << std::endl;
  #endif
  #ifdef KOKKOS_ENABLE_CUDA
    idfx::cout << "Input: Kokkos CUDA target ENABLED." << std::endl;
  #endif
  #ifdef KOKKOS_ENABLE_OPENMP
    idfx::cout << "Input: Kokkos OpenMP ENABLED." << std::endl;
  #endif
}

// This routine is called whenever a specific OS signal is caught
void Input::signalHandler(int signum) {
  idfx::cout << std::endl << "Input: Caught interrupt " << signum << std::endl;
  abortRequested=true;
}

void Input::CheckForStopFile() {
  // Check whether a file "stop" has been created in directory. If so, raise the abort flag
  // Look for stop file every 5 seconds to avoid overloading the file system
  if(idfx::prank==0) {
    std::string filename = std::string("stop");
    if(timer.seconds()-lastStopFileCheck>5.0) {
      lastStopFileCheck = timer.seconds();
      std::ifstream f(filename);
      if(f.good()) {
        // File exists, delete it and raise the flag
        std::remove(filename.c_str());
        abortRequested = true;
        idfx::cout << std::endl << "Input: Caught stop file command" << std::endl;
      }
    }
  }
}

bool Input::CheckForAbort() {
  idfx::pushRegion("Input::CheckForAbort");
  // Check whether an abort has been requesested
  // When MPI is present, we abort whenever one process got the signal
  CheckForStopFile();
#ifdef WITH_MPI
  int abortValue{0};
  bool returnValue{false};
  if(abortRequested) abortValue = 1;

  MPI_Bcast(&abortValue, 1, MPI_INT, 0, MPI_COMM_WORLD);
  returnValue = abortValue > 0;
  if(returnValue) idfx::cout << "Input: CheckForAbort: abort has been requested." << std::endl;
  idfx::popRegion();
  return(returnValue);
#else
  if(abortRequested) idfx::cout << "Input: CheckForAbort: abort has been requested." << std::endl;
  idfx::popRegion();
  return(abortRequested);
#endif
}

// Get a string in a block, parameter, position of the file
std::string Input::GetString(std::string blockName, std::string paramName, int num) {
  IDEFIX_DEPRECATED("Input::GetString is deprecated. Use Input::Get<std::string> instead");
  return(Get<std::string>(blockName, paramName, num));
}

// Get a real number in a block, parameter, position of the file
real Input::GetReal(std::string blockName, std::string paramName, int num) {
  IDEFIX_DEPRECATED("Input::GetReal is deprecated. Use Input::Get<real> instead");
  return(Get<real>(blockName, paramName, num));
}

// Get an integer number in a block, parameter, position of the file
int Input::GetInt(std::string blockName, std::string paramName, int num) {
  IDEFIX_DEPRECATED("Input::GetInt is deprecated. Use Input::Get<int> instead");
  return(Get<int>(blockName, paramName, num));
}

// Check that an entry is present in the ini file.
// If yes, return the number of parameters for given entry
int Input::CheckEntry(std::string blockName, std::string paramName) {
  int result=-1;
  IdefixInputContainer::iterator block = inputParameters.find(blockName);
  if(block != inputParameters.end()) {
    // Block exists
    IdefixBlockContainer::iterator param = block->second.find(paramName);
    if(param != block->second.end()) {
      // Parameter exist
      result = param->second.size();
    }
  }
  return(result);
}

// Check that a block is present in the ini file.
// If yes, return true
bool Input::CheckBlock(std::string blockName) {
  bool result = false;
  IdefixInputContainer::iterator block = inputParameters.find(blockName);
  if(block != inputParameters.end()) {
    // Block exists
    result = true;
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
  PrintVersion();
  idfx::cout << std::endl;
  idfx::cout << std::endl;
}

void Input::PrintOptions() {
  idfx::cout << "List of valid arguments:" << std::endl << std::endl;
  #ifdef WITH_MPI
    idfx::cout << " -dec "
    D_SELECT(<< " nx1 ",
             << " nx1 nx2",
             << " nx1 nx2 nx3")
            << std::endl;
    idfx::cout << "         Force an mpi domain decomposition with n processes in each direction"
               << std::endl;
  #endif
  idfx::cout << " -restart n" << std::endl;
  idfx::cout << "         Restart from dumpfile n. If n is ommited, Idefix restart from the latest"
             << " generated dump file." << std::endl;
  idfx::cout << " -i xxx" << std::endl;
  idfx::cout << "         Use the input file xxx instead of the default idefix.ini" << std::endl;
  idfx::cout << " -maxcycles n" << std::endl;
  idfx::cout << "         Perform at most n integration cycles." << std::endl;
  idfx::cout << " -force_init" << std::endl;
  idfx::cout << "         Call initial conditions before reading dump file ";
  idfx::cout << "(this has no effect if -restart is not also passed)" << std::endl;
  idfx::cout << " -nowrite" << std::endl;
  idfx::cout << "         Do not generate any output file." << std::endl;
  idfx::cout << " -nolog" << std::endl;
  idfx::cout << "         Do not write any log file." << std::endl;
  idfx::cout << " -profile" << std::endl;
  idfx::cout << "         Enable on-the-fly performance profiling." << std::endl;
  idfx::cout << " -Werror" << std::endl;
  idfx::cout << "         Consider warnings as errors." << std::endl;
  idfx::cout << " -v/-version" << std::endl;
  idfx::cout << "         Show Idefix and kokkos version." << std::endl;
  idfx::cout << " -h/-help" << std::endl;
  idfx::cout << "         Show this message." << std::endl;
}

void Input::PrintVersion() {
  idfx::cout << "              Idefix version " << IDEFIX_VERSION << std::endl;
  idfx::cout << "              Built against Kokkos " << KOKKOS_VERSION << std::endl;
  idfx::cout << "              Compiled on " << __DATE__ <<  " at " << __TIME__ << std::endl;
}
