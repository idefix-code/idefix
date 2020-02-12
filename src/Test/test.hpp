#ifndef TEST_HPP
#define TEST_HPP
#include "idefix.hpp"


class Test {
public:

    // Constructor
    Test(DataBlock&);

    // Default constructor
    Test();

    // Make the test
    void MakeTest(int, int );

private:
    DataBlock data;
    DataBlockHost dataHost;
};

#endif