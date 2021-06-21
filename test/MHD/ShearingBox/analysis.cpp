#include "analysis.hpp"
#include "idefix.hpp"
#include <iostream>
#include <fstream>

Analysis::Analysis(Input &input, Grid &grid, DataBlock &data, Output &output, std::string filename) {
      this->d = new DataBlockHost(data);
      this->grid = &grid;
      this->filename = filename;
      this->shear = data.hydro.sbS;
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
    file << std::setw(col_width) << "rho";
    file << std::setw(col_width) << "prs";
    file << std::setw(col_width) << "rhovx2";
    file << std::setw(col_width) << "rhovy2";
    file << std::setw(col_width) << "rhovz2";
    file << std::setw(col_width) << "rhovxvy";
    file << std::setw(col_width) << "Bx";
    file << std::setw(col_width) << "By";
    file << std::setw(col_width) << "Bz";
    file << std::setw(col_width) << "BxBy";
    file << std::endl;
    file.close();
  }
}

void Analysis::PerformAnalysis(DataBlock &data) {
  idfx::pushRegion("Analysis::PerformAnalysis");
  d->SyncFromDevice();
  int fields[3];
  int nfields;
  if(idfx::prank==0) {
    file.open(filename, std::ios::app);
    file.precision(precision);
  }

  WriteField(data.t);

  fields[0] = RHO;
  WriteField(Average(1, fields));

  fields[0] = PRS;
  WriteField(Average(1, fields));

  fields[0] = RHO;
  fields[1] = VX1;
  fields[2] = VX1;
  WriteField(Average(3, fields));

  fields[0] = RHO;
  fields[1] = VX2;
  fields[2] = VX2;
  WriteField(Average(3, fields));

  fields[0] = RHO;
  fields[1] = VX3;
  fields[2] = VX3;
  WriteField(Average(3, fields));

  fields[0] = RHO;
  fields[1] = VX1;
  fields[2] = VX2;
  WriteField(Average(3, fields));

  fields[0] = BX1;
  WriteField(Average(1, fields));

  fields[0] = BX2;
  WriteField(Average(1, fields));

  fields[0] = BX3;
  WriteField(Average(1, fields));

  fields[0] = BX1;
  fields[1] = BX2;
  WriteField(Average(2, fields));

  if(idfx::prank==0) {
    file << std::endl;
    file.close();
  }
  idfx::popRegion();
}
