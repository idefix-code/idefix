#ifndef GLOBAL_HPP
#define GLOBAL_HPP

namespace idfx {
    int initialize();   // Initialisation routine for idefix
    class IdefixOstream;

    extern int prank;       // parallel rank
    extern int psize;
    extern IdefixOstream cout;  // custom cout for idefix

    void pushRegion(const std::string&);
    void popRegion();


}

class idfx::IdefixOstream
{
public:
 
  void init(int);
  // for regular output of variables and stuff
  template<typename T> IdefixOstream& operator<<(const T& something)
  {
    if(toscreen) std::cout << something;
    my_fstream << something;
    return *this;
  }
  // for manipulators like std::endl
  typedef std::ostream& (*stream_function)(std::ostream&);
  IdefixOstream& operator<<(stream_function func)
  {
    if(toscreen) func(std::cout);
    func(my_fstream);
    return *this;
  }
private:
  std::ofstream my_fstream;
  bool toscreen;
};

#endif