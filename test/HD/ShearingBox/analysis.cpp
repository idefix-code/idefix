#include "analysis.hpp"
#include "idefix.hpp"
#include "fluid.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

Analysis::Analysis(Input &input, Grid &grid, DataBlock &data, Output &output, std::string filename) {
      this->d = new DataBlockHost(data);
      this->grid = &grid;
      this->filename = filename;
      this->shear = data.hydro->sbS;
      this->precision = 10;
}

/* **************************************************************** */
double Analysis::Average(const int nfields, int fields[])
    /*
     * compute the weighted average: int dphi dz rho *infield/int dphi dz rho
     *
     **************************************************************** */
{
  real outfield = 0;

  for(int k = d->beg[KDIR]; k < d->end[KDIR] ; k++) {
    for(int j = d->beg[JDIR]; j < d->end[JDIR] ; j++) {
      for(int i = d->beg[IDIR]; i < d->end[IDIR] ; i++) {
        real q=1.0;
        for(int n=0 ; n < nfields ; n++) {
          // Substract Keplerian flow if vphi
          if(fields[n]==VX2) {
            q = q*(d->Vc(fields[n],k,j,i)-d->x[IDIR](i)*this->shear);
          }
          else{
            q = q*d->Vc(fields[n],k,j,i);
          }
        }
        outfield += q;
      }
    }
  }

    // Reduce
#ifdef WITH_MPI
  real reducedValue;
  MPI_Reduce(&outfield, &reducedValue, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  outfield = reducedValue;
#endif

  outfield = outfield / ((double) grid->np_int[IDIR] * grid->np_int[JDIR] * grid->np_int[KDIR]);


  return outfield;
}

/* **************************************************************** */
double Analysis::ShwaveAmplitude(const int field,
                                 const int nx,
                                 const int ny,
                                 const int nz,
                                 const real t)
    /*
     * compute the weighted average: int dphi dz rho *infield/int dphi dz rho
     *
     **************************************************************** */
{

  real q =0.0;
  for(int k = d->beg[KDIR]; k < d->end[KDIR] ; k++) {
    for(int j = d->beg[JDIR]; j < d->end[JDIR] ; j++) {
      for(int i = d->beg[IDIR]; i < d->end[IDIR] ; i++) {
        real x = d->x[IDIR](i);
        real y = d->x[JDIR](j);
        real z = d->x[KDIR](k);
        real wave = sin(2.0*M_PI*(  (nx-ny*shear*t)*x
                                  + ny*y
                                  + nz*z));
        q += wave*d->Vc(field,k,j,i);
      }
    }
  }

    // Reduce
#ifdef WITH_MPI
  real reducedValue;
  MPI_Reduce(&q, &reducedValue, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  q = reducedValue;
#endif

  q = q / ((double) grid->np_int[IDIR] * grid->np_int[JDIR] * grid->np_int[KDIR]);


  return q;
}

/* **************************************************************** */
void Analysis::WriteField(double data) {
/*
 * Write a global profile to a file
 *
 *
 **************************************************************** */
  if(idfx::prank==0) {
    int col_width = precision + 10;
    this->file << std::scientific << std::setw(col_width) << data;
  }
  return ;
}


void Analysis::ResetAnalysis() {
  GridHost gh(*this->grid);
  gh.SyncFromDevice();
  int col_width = precision + 10;
  if(idfx::prank==0) {
    file.open(filename, std::ios::trunc);
    file << std::setw(col_width) << "t";
    file << std::setw(col_width) << "vx2";
    file << std::setw(col_width) << "vy2";
    file << std::setw(col_width) << "vz2";
    file << std::setw(col_width) << "vx";
    file << std::setw(col_width) << "vy";
    file << std::setw(col_width) << "vz";
    file << std::endl;
    file.close();
  }
}

void Analysis::PerformAnalysis(DataBlock &data) {
  idfx::pushRegion("Analysis::PerformAnalysis");
  d->SyncFromDevice();
  int fields[3];
  if(idfx::prank==0) {
    file.open(filename, std::ios::app);
    file.precision(precision);
  }

  const int nx0 =  data.beg[IDIR] + data.np_int[IDIR]/2;
  const int ny0 =  data.beg[JDIR] + 3*data.np_int[JDIR]/4;
  const int nz0 = data.beg[KDIR];

  WriteField(data.t);

  fields[0] = VX1;
  fields[1] = VX1;
  WriteField(Average(2, fields));

  fields[0] = VX2;
  fields[1] = VX2;
  WriteField(Average(2, fields));

  fields[0] = VX3;
  fields[1] = VX3;
  WriteField(Average(2, fields));

  WriteField(ShwaveAmplitude(VX1, 0, 1, 4, data.t) );
  WriteField(ShwaveAmplitude(VX2, 0, 1, 4, data.t) );
  WriteField(ShwaveAmplitude(VX3, 0, 1, 4, data.t) );


  if(idfx::prank==0) {
    file << std::endl;
    file.close();
  }
  idfx::popRegion();
}
