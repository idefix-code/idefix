// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include <vector>
#include <string>

#include "column.hpp"
#include "idefix.hpp"
#include "input.hpp"
#include "output.hpp"
#include "grid.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"

Column::Column(int dir, int sign, DataBlock *data)
                : direction(dir), sign(sign) {
  idfx::pushRegion("Column::Column");
  this->np_tot = data->np_tot;
  this->np_int = data->np_int;
  this->beg = data->beg;

  if(dir>= DIMENSIONS || dir < IDIR) IDEFIX_ERROR("Unknown direction for Column constructor");

  // Allocate the array on which we do the average
  this->ColumnArray = IdefixArray3D<real>("ColumnArray",np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);

  // Elementary volumes and area
  this->Volume = data->dV;
  this->Area = data->A[dir];

  // allocate helper  array
  if(dir == IDIR) {
    localSum = IdefixArray2D<real>("localSum",np_tot[KDIR], np_tot[JDIR]);
  }
  if(dir == JDIR) {
    localSum = IdefixArray2D<real>("localSum",np_tot[KDIR], np_tot[IDIR]);
  }
  if(dir == KDIR) {
    localSum = IdefixArray2D<real>("localSum",np_tot[JDIR], np_tot[IDIR]);
  }
  #ifdef WITH_MPI
  // Create sub-MPI communicator dedicated to scan
    int remainDims[3] = {false, false, false};
    remainDims[dir] = true;
    MPI_Cart_sub(data->mygrid->CartComm, remainDims, &this->ColumnComm);
    MPI_Comm_rank(this->ColumnComm, &this->MPIrank);
    MPI_Comm_size(this->ColumnComm, &this->MPIsize);

    // create MPI class for boundary Xchanges
    int ntarget = 0;
    std::vector<int> mapVars;
    mapVars.push_back(ntarget);

    this->mpi.Init(data->mygrid, mapVars, data->nghost.data(), data->np_int.data());
    this->nproc = data->mygrid->nproc;
  #endif
  idfx::popRegion();
}

void Column::ComputeColumn(IdefixArray4D<real> in, const int var) {
  idfx::pushRegion("Column::ComputeColumn");
  const int nk = np_int[KDIR];
  const int nj = np_int[JDIR];
  const int ni = np_int[IDIR];
  const int kb = beg[KDIR];
  const int jb = beg[JDIR];
  const int ib = beg[IDIR];
  const int ke = kb+nk;
  const int je = jb+nj;
  const int ie = ib+ni;

  const int direction = this->direction;
  auto column = this->ColumnArray;
  auto dV = this->Volume;
  auto A = this->Area;
  auto localSum = this->localSum;

  if(direction==IDIR) {
    // Inspired from loop.hpp
    Kokkos::parallel_for("ColumnX1", team_policy (nk*nj, Kokkos::AUTO),
      KOKKOS_LAMBDA (member_type team_member) {
        int k = team_member.league_rank() / nj;
        int j = team_member.league_rank() - k*nj + jb;
        k += kb;
        Kokkos::parallel_scan(Kokkos::TeamThreadRange<>(team_member,ib,ie),
          [=] (int i, real &partial_sum, bool is_final) {
            partial_sum += in(var,k,j,i)*dV(k,j,i) / (0.5*(A(k,j,i)+A(k,j,i+1)));
            if(is_final) column(k,j,i) = partial_sum;
            //partial_sum += in(var,k,j,i)*dV(k,j,i) / (0.5*(A(k,j,i)+A(k,j,i+1)));
          });
      });
    }
    if(direction==JDIR) {
      // Inspired from loop.hpp
      Kokkos::parallel_for("ColumnX2", team_policy (nk*ni, Kokkos::AUTO),
        KOKKOS_LAMBDA (member_type team_member) {
          int k = team_member.league_rank() / ni;
          int i = team_member.league_rank() - k*ni + ib;
          k += kb;
          Kokkos::parallel_scan(Kokkos::TeamThreadRange<>(team_member,jb,je),
            [=] (int j, real &partial_sum, bool is_final) {
              partial_sum += in(var,k,j,i)*dV(k,j,i) / (0.5*(A(k,j,i)+A(k,j+1,i)));
              if(is_final) column(k,j,i) = partial_sum;
          });
      });
    }
    if(direction==KDIR) {
      // Inspired from loop.hpp
      Kokkos::parallel_for("ColumnX2", team_policy (nj*ni, Kokkos::AUTO),
        KOKKOS_LAMBDA (member_type team_member) {
          int j = team_member.league_rank() / ni;
          int i = team_member.league_rank() - j*ni + ib;
          j += jb;
          Kokkos::parallel_scan(Kokkos::TeamThreadRange<>(team_member,kb,ke),
            [=] (int k, real &partial_sum, bool is_final) {
              partial_sum += in(var,k,j,i)*dV(k,j,i) / (0.5*(A(k,j,i)+A(k+1,j,i)));
              if(is_final) column(k,j,i) = partial_sum;
          });
      });
    }

    #ifdef WITH_MPI
    // Load the current sum
      int dst,src;
      MPI_Cart_shift(this->ColumnComm,0, 1, &src, &dst);
      int size = localSum.extent(0)*localSum.extent(1);
      if(MPIrank>0) {
        MPI_Status status;
        // Get the cumulative sum from previous processes
        Kokkos::fence();
        MPI_Recv(localSum.data(), size, realMPI, src, 20, ColumnComm, &status);
        // Add this to our cumulative sum
        idefix_for("Addsum",kb,ke,jb,je,ib,ie,
          KOKKOS_LAMBDA(int k, int j, int i) {
            if(direction == IDIR) column(k,j,i) += localSum(k,j);
            if(direction == JDIR) column(k,j,i) += localSum(k,i);
            if(direction == KDIR) column(k,j,i) += localSum(j,i);
        });
      } // rank =0
      // Get the final sum
      if(MPIrank < MPIsize-1) {
        if(direction==IDIR) {
          idefix_for("Loadsum",kb,ke,jb,je,
            KOKKOS_LAMBDA(int k, int j) {
              localSum(k,j) = column(k,j,ie-1);
          });
        }
        if(direction==JDIR) {
          idefix_for("Loadsum",kb,ke,ib,ie,
            KOKKOS_LAMBDA(int k, int i) {
              localSum(k,i) = column(k,je-1,i);
          });
        }
        if(direction==KDIR) {
          idefix_for("Loadsum",jb,je,ib,ie,
            KOKKOS_LAMBDA(int j, int i) {
              localSum(j,i) = column(ke-1,j,i);
          });
        }
        // And send it
        Kokkos::fence();
        MPI_Send(localSum.data(),size, realMPI, dst, 20, ColumnComm);
      } // MPIrank small enough
    #endif
    // If we need it backwards
    if(sign<0) {
      // Compute total cumulative sum in the last subdomain
      if(direction == IDIR) {
        idefix_for("Loadsum",kb,ke,jb,je,
          KOKKOS_LAMBDA(int k, int j) {
            localSum(k,j) = column(k,j,ie-1) + in(var,k,j,ie-1) * dV(k,j,ie-1)
                            / (0.5*(A(k,j,ie-1)+A(k,j,ie)));
        });
      }
      if(direction == JDIR) {
        idefix_for("Loadsum",kb,ke,ib,ie,
          KOKKOS_LAMBDA(int k, int i) {
            localSum(k,i) = column(k,je-1,i) + in(var,k,je-1,i) * dV(k,je-1,i)
                            / (0.5*(A(k,je-1,i)+A(k,je,i)));
        });
      }
      if(direction == KDIR) {
        idefix_for("Loadsum",jb,je,ib,ie,
          KOKKOS_LAMBDA(int j, int i) {
            localSum(j,i) = column(ke-1,j,i) + in(var,ke-1,j,i) * dV(ke-1,j,i)
                            / (0.5*(A(ke-1,j,i)+A(ke,j,i)));
        });
      }

      #ifdef WITH_MPI
      Kokkos::fence();
      MPI_Bcast(localSum.data(),size, realMPI, MPIsize-1, ColumnComm);
      #endif
      // All substract the local column from the full column

      idefix_for("Subsum",kb,ke,jb,je,ib,ie,
        KOKKOS_LAMBDA(int k, int j, int i) {
          if(direction==IDIR) column(k,j,i) = localSum(k,j)-column(k,j,i);
          if(direction==JDIR) column(k,j,i) = localSum(k,i)-column(k,j,i);
          if(direction==KDIR) column(k,j,i) = localSum(j,i)-column(k,j,i);
      });
    } // sign <0
    // Xchange boundary elements when using MPI to ensure that column
    // density in the ghost zones are coherent
    #ifdef WITH_MPI
      // Create a 4D array that contains our column data
      IdefixArray4D<real> arr4D(column.data(), 1, this->np_tot[KDIR],
                                                        this->np_tot[JDIR],
                                                        this->np_tot[IDIR]);

      for(int dir = 0 ; dir < DIMENSIONS ; dir++) {
        // MPI Exchange data when needed
        if(nproc[dir]>1) {
          switch(dir) {
            case 0:
              this->mpi.ExchangeX1(arr4D);
              break;
            case 1:
              this->mpi.ExchangeX2(arr4D);
              break;
            case 2:
              this->mpi.ExchangeX3(arr4D);
              break;
          }
        }
      }
    #endif
  idfx::popRegion();
}

void Column::ComputeColumn(IdefixArray3D<real> in) {
  // 4D alias
  IdefixArray4D<real> arr4D(in.data(), 1, in.extent(2), in.extent(1), in.extent(0));
  return this->ComputeColumn(arr4D,0);
}
