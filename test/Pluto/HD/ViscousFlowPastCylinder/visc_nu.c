/* /////////////////////////////////////////////////////////////////// */
/*! \file
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stdio.h>
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
             double *nu1, double *nu2)
/*!
 * \brief Calculate first and second viscosity coefficients as functions of data and coordinates
 *
 *    \param [in]      v  pointer to data array containing cell-centered quantities
 *    \param [in]      x1 real, coordinate value
 *    \param [in]      x2 real, coordinate value
 *    \param [in]      x3 real, coordinate value
 *    \param [in, out] nu1 pointer to first viscous coefficient
 *    \param [in, out] nu2 pointer to second viscous coefficient
 *    \return This function has no return value.
 * ************************************************************************** */

{
 *nu1 = g_inputParam[NU_VISC];
 *nu2 = 0.0;
}
