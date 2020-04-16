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
        std::cout << "------------------------------------------------------------------------------" << std::endl;
        std::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        std::cout << ErrorMessage.str() << std::endl;
        std::cout << "------------------------------------------------------------------------------" << std::endl;
    }
    else
    {
        std::cout << "------------------------------------------------------------------------------" << std::endl;
        std::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        std::cout << ErrorMessage.str() << std::endl;
        std::cout << "------------------------------------------------------------------------------" << std::endl;
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
        std::cout << "------------------------------------------------------------------------------" << std::endl;
        std::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        std::cout << ErrorMessage << std::endl;
        std::cout << "------------------------------------------------------------------------------" << std::endl;
    }
    else
    {
        std::cout << "------------------------------------------------------------------------------" << std::endl;
        std::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile << " line " << ErrorLine << std::endl;
        std::cout << ErrorMessage << std::endl;
        std::cout << "------------------------------------------------------------------------------" << std::endl;
        exit(1);
    }
}
