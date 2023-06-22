// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line "
               << ErrorLine << std::endl;
    idfx::cerr << ErrorMessage.str() << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    if(idfx::warningsAreErrors) {
      idfx::cerr << "Warnings are considered as errors" << std::endl;
      idfx::safeExit(1);
    }
  } else if (ErrorType == ERROR_DEPRECATED) {
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "## DEPRECATED call in function " << ErrorFunction << " file " << ErrorFile
               <<  std::endl
               << "## This function will be removed in the next Idefix release." << std::endl;
    idfx::cerr << ErrorMessage.str() << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    if(idfx::warningsAreErrors) {
      idfx::cerr << "Warnings are considered as errors" << std::endl;
      idfx::safeExit(1);
    }
  } else {
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile
               << " line " << ErrorLine << std::endl;
    idfx::cerr << ErrorMessage.str() << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::safeExit(1);
  }
}

void ErrorHandler(const int ErrorType,
                  std::string ErrorMessage,
                  std::string ErrorFunction,
                  const int ErrorLine,
                  std::string ErrorFile) {
  if (ErrorType == ERROR_WARNING) {
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "## WARNING in function " << ErrorFunction << " file " << ErrorFile << " line "
               << ErrorLine << std::endl;
    idfx::cerr << ErrorMessage << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    if(idfx::warningsAreErrors) {
      idfx::cerr << "Warnings are considered as errors" << std::endl;
      idfx::safeExit(1);
    }
  } else if (ErrorType == ERROR_DEPRECATED) {
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "## DEPRECATED call in function " << ErrorFunction << " file " << ErrorFile
               <<  std::endl
               << "## This function will be removed in the next Idefix release." << std::endl;
    idfx::cerr << ErrorMessage << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    if(idfx::warningsAreErrors) {
      idfx::cerr << "Warnings are considered as errors" << std::endl;
      idfx::safeExit(1);
    }
  } else {
    idfx::cerr << std::endl
               << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::cerr << "### FATAL ERROR in function " << ErrorFunction << " file " << ErrorFile
               << " line " << ErrorLine << std::endl;
    idfx::cerr << ErrorMessage << std::endl;
    idfx::cerr << "------------------------------------------------------------------------------"
               << std::endl;
    idfx::safeExit(1);
  }
}
