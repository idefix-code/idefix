// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_
#include <string>

namespace idfx {
int initialize();   // Initialisation routine for idefix
class IdefixOstream;
class Profiler;

extern int prank;     // parallel rank
extern int psize;
extern IdefixOstream cout;  // custom cout for idefix
extern Profiler prof;     // profiler (for memory usage)
extern double mpiCallsTimer;            //< time significant MPI calls
extern LoopPattern defaultLoopPattern;  //< default loop patterns (for idefix_for loops)

void pushRegion(const std::string&);
void popRegion();
} // namespace idfx

class idfx::IdefixOstream {
 public:
  void init(int);
  // for regular output of variables and stuff
  template<typename T> IdefixOstream& operator<<(const T& something) {
    if(toscreen) std::cout << something;
    my_fstream << something;
    return *this;
  }
  // for manipulators like std::endl
  typedef std::ostream& (*stream_function)(std::ostream&);
  IdefixOstream& operator<<(stream_function func) {
    if(toscreen) func(std::cout);
    func(my_fstream);
    return *this;
  }
 private:
  std::ofstream my_fstream;
  bool toscreen;
};

#endif // GLOBAL_HPP_
