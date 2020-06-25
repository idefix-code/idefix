#include <cstdio>
#include "idefix.hpp"
#include "error.hpp"

void ErrorHandler(const int ErrorType,
                  std::stringstream &ErrorMessage,
                  std::string ErrorFunction,
                  const int ErrorLine,
                  std::string ErrorFile)
{

    if (ErrorType == ERROR_WARNING)
    {
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        idfx::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        idfx::cout << ErrorMessage.str() << std::endl;
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
    }
    else
    {
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        idfx::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        idfx::cout << ErrorMessage.str() << std::endl;
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        #ifdef WITH_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
        #endif
        exit(1);
    }
}

void ErrorHandler(const int ErrorType,
                  std::string ErrorMessage,
                  std::string ErrorFunction,
                  const int ErrorLine,
                  std::string ErrorFile)
{

    if (ErrorType == ERROR_WARNING)
    {
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        idfx::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        idfx::cout << ErrorMessage << std::endl;
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
    }
    else
    {
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        idfx::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        idfx::cout << ErrorMessage << std::endl;
        idfx::cout << "------------------------------------------------------------------------------" << std::endl;
        #ifdef WITH_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
        #endif
        exit(1);
    }
}
