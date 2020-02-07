#ifndef TEST_HPP
#define TEST_HPP
#include "idefix.hpp"


class Test {
public:

    // Constructor
    Test(Data&);

    // Default constructor
    Test();

    // Make the test
    void MakeTest(Grid &, int, int );

private:
    Data data;
    DataHost dataHost;
};

#endif