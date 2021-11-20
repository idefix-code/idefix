// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <cstdio>
#include <string>
#include "idefix.hpp"
#include "error.hpp"

void ErrorHandler(const int ErrorType,
                  std::stringstream &ErrorMessage,
                  std::string ErrorFunction,
                  const int ErrorLine,
                  std::string ErrorFile) {
  if (ErrorType == ERROR_WARNING) {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line "
               << ErrorLine << std::endl;
    idfx::cout << ErrorMessage.str() << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
  } else if (ErrorType == ERROR_DEPRECATED) {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "## DEPRECATED call in function " << ErrorFunction << " file " << ErrorFile
               <<  std::endl
               << "## This function will be removed in the next Idefix release." << std::endl;
    idfx::cout << ErrorMessage.str() << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
  } else {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile
               << " line " << ErrorLine << std::endl;
    idfx::cout << ErrorMessage.str() << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
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
                  std::string ErrorFile) {
  if (ErrorType == ERROR_WARNING) {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line "
               << ErrorLine << std::endl;
    idfx::cout << ErrorMessage << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
  } else if (ErrorType == ERROR_DEPRECATED) {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "## DEPRECATED call in function " << ErrorFunction << " file " << ErrorFile
               <<  std::endl
               << "## This function will be removed in the next Idefix release." << std::endl;
    idfx::cout << ErrorMessage << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
  } else {
    idfx::cout << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cout << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile
               << " line " << ErrorLine << std::endl;
    idfx::cout << ErrorMessage << std::endl;
    idfx::cout << "------------------------------------------------------------------------------"
               << std::endl;
    #ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
    #endif
    exit(1);
  }
}
