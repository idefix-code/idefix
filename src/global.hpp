// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_
#include <string>

namespace idfx {
int initialize();   // Initialisation routine for idefix
class IdefixOstream;

extern int prank;     // parallel rank
extern int psize;
extern IdefixOstream cout;  // custom cout for idefix
extern double mpiTimer;   // Measure MPI perfs as you go

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
