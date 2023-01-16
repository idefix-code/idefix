// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef ERROR_HPP_
#define ERROR_HPP_

#include <string>


#define ERROR_WARNING 1
#define ERROR_CRITICAL 2
#define ERROR_DEPRECATED 3


#define IDEFIX_ERROR(ERROR_MESSAGE) ErrorHandler(ERROR_CRITICAL, ERROR_MESSAGE, \
                                         std::string(__func__) , __LINE__, std::string(__FILE__))
#define IDEFIX_WARNING(WARNING_MESSAGE) ErrorHandler(ERROR_WARNING, WARNING_MESSAGE, \
                                         std::string(__func__) , __LINE__, std::string(__FILE__))
#define IDEFIX_DEPRECATED(WARNING_MESSAGE) ErrorHandler(ERROR_DEPRECATED, WARNING_MESSAGE, \
                                         std::string(__func__) , __LINE__, std::string(__FILE__))

void ErrorHandler (const int, std::stringstream &, std::string , const int, std::string );

void ErrorHandler (const int, std::string, std::string, const int, std::string );

#endif // ERROR_HPP_
