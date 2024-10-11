// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_
#include <iostream>
#include <string>
#include <vector>
#include "arrays.hpp"

namespace idfx {
int initialize();   // Initialisation routine for idefix
real randm();      // Custom random number generator
void safeExit(int );       // Exit the code
class IdefixOutStream;
class IdefixErrStream;
class Profiler;

extern int prank;                       //< parallel rank
extern int psize;
extern IdefixOutStream cout;              //< custom cout for idefix
extern IdefixErrStream cerr;              //< custom cerr for idefix
extern Profiler prof;                   //< profiler (for memory & performance usage)
extern double mpiCallsTimer;            //< time significant MPI calls
extern LoopPattern defaultLoopPattern;  //< default loop patterns (for idefix_for loops)
extern bool warningsAreErrors;    //< whether warnings should be considered as errors

void pushRegion(const std::string&);
void popRegion();

template<typename T>
IdefixArray1D<T> ConvertVectorToIdefixArray(std::vector<T> &inputVector) {
  IdefixArray1D<T> outArr = IdefixArray1D<T>("Vector",inputVector.size());
  IdefixHostArray1D<T> outArrHost;
  outArrHost = Kokkos::create_mirror_view(outArr);
  for(int i = 0; i < inputVector.size() ; i++) {
    outArrHost(i) = inputVector[i];
  }
  Kokkos::deep_copy(outArr, outArrHost);
  return(outArr);
}

} // namespace idfx

class idfx::IdefixOutStream {
 public:
  void init(int);
  void enableLogFile();
  // for regular output of variables and stuff
  template<typename T> IdefixOutStream& operator<<(const T& something) {
    if(toscreen) std::cout << something;
    if(logFileEnabled) my_fstream << something;
    return *this;
  }
  // for manipulators like std::endl
  typedef std::ostream& (*stream_function)(std::ostream&);
  IdefixOutStream& operator<<(stream_function func) {
    if(toscreen) func(std::cout);
    if(logFileEnabled) func(my_fstream);
    return *this;
  }
 private:
  std::ofstream my_fstream;
  bool toscreen;
  bool logFileEnabled{false};   //< whether streams are also written to a log file
};

class idfx::IdefixErrStream {
 public:
  // for error output of variables and stuff
  template<typename T> IdefixErrStream& operator<<(const T& something) {
    std::cerr << something;
    return *this;
  }
  // for manipulators like std::endl
  typedef std::ostream& (*stream_function)(std::ostream&);
  IdefixErrStream& operator<<(stream_function func) {
    func(std::cerr);
    return *this;
  }
};

#endif // GLOBAL_HPP_
