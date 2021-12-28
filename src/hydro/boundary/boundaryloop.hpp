// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#ifndef HYDRO_BOUNDARY_BOUNDARYLOOP_HPP_
#define HYDRO_BOUNDARY_BOUNDARYLOOP_HPP_

#include <string>
#include "idefix.hpp"

template <typename Function>
inline void HydroBoundary::BoundaryForAll(
  const std::string & name,
  const int &dir,
  const BoundarySide &side,
  Function function) {
    const int nxi = data->np_int[IDIR];
    const int nxj = data->np_int[JDIR];
    const int nxk = data->np_int[KDIR];

    const int ighost = data->nghost[IDIR];
    const int jghost = data->nghost[JDIR];
    const int kghost = data->nghost[KDIR];

    // Boundaries of the loop
    const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
    const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR];
    const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
    const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR];
    const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
    const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR];

    idefix_for(name, 0, NVAR, kbeg, kend, jbeg, jend, ibeg, iend, function);
}

template <typename Function>
inline void HydroBoundary::BoundaryFor(
  const std::string & name,
  const int &dir,
  const BoundarySide &side,
  Function function) {
    const int nxi = data->np_int[IDIR];
    const int nxj = data->np_int[JDIR];
    const int nxk = data->np_int[KDIR];

    const int ighost = data->nghost[IDIR];
    const int jghost = data->nghost[JDIR];
    const int kghost = data->nghost[KDIR];

    // Boundaries of the loop
    const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
    const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR];
    const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
    const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR];
    const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
    const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR];



    idefix_for(name, kbeg, kend, jbeg, jend, ibeg, iend, function);
}

template <typename Function>
inline void HydroBoundary::BoundaryForX1s(
  const std::string & name,
  const int &dir,
  const BoundarySide &side,
  Function function) {
    const int nxi = data->np_int[IDIR]+1;
    const int nxj = data->np_int[JDIR];
    const int nxk = data->np_int[KDIR];

    const int ighost = data->nghost[IDIR];
    const int jghost = data->nghost[JDIR];
    const int kghost = data->nghost[KDIR];

    // Boundaries of the loop
    const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
    const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR]+1;
    const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
    const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR];
    const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
    const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR];

    idefix_for(name, kbeg, kend, jbeg, jend, ibeg, iend, function);
}

template <typename Function>
inline void HydroBoundary::BoundaryForX2s(
  const std::string & name,
  const int &dir,
  const BoundarySide &side,
  Function function) {
    const int nxi = data->np_int[IDIR];
    const int nxj = data->np_int[JDIR]+1;
    const int nxk = data->np_int[KDIR];

    const int ighost = data->nghost[IDIR];
    const int jghost = data->nghost[JDIR];
    const int kghost = data->nghost[KDIR];

    // Boundaries of the loop
    const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
    const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR];
    const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
    const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR]+1;
    const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
    const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR];

    idefix_for(name, kbeg, kend, jbeg, jend, ibeg, iend, function);
}

template <typename Function>
inline void HydroBoundary::BoundaryForX3s(
  const std::string & name,
  const int &dir,
  const BoundarySide &side,
  Function function) {
    const int nxi = data->np_int[IDIR];
    const int nxj = data->np_int[JDIR];
    const int nxk = data->np_int[KDIR]+1;

    const int ighost = data->nghost[IDIR];
    const int jghost = data->nghost[JDIR];
    const int kghost = data->nghost[KDIR];

    // Boundaries of the loop
    const int ibeg = (dir == IDIR) ? side*(ighost+nxi) : 0;
    const int iend = (dir == IDIR) ? ighost + side*(ighost+nxi) : data->np_tot[IDIR];
    const int jbeg = (dir == JDIR) ? side*(jghost+nxj) : 0;
    const int jend = (dir == JDIR) ? jghost + side*(jghost+nxj) : data->np_tot[JDIR];
    const int kbeg = (dir == KDIR) ? side*(kghost+nxk) : 0;
    const int kend = (dir == KDIR) ? kghost + side*(kghost+nxk) : data->np_tot[KDIR]+1;

    idefix_for(name, kbeg, kend, jbeg, jend, ibeg, iend, function);
}




#endif // HYDRO_BOUNDARY_BOUNDARYLOOP_HPP_
