// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef UTILS_LOOKUPTABLE_HPP_
#define UTILS_LOOKUPTABLE_HPP_

#include <string>
#include <vector>
#include "idefix.hpp"
#include "lookupTable.hpp"
#include "npy.hpp"

template <const int kDim>
class LookupTable {
 public:
  LookupTable() = default;
  LookupTable(std::string filename, char delimiter, bool errorIfOutOfBound = true);
  LookupTable(std::vector<std::string> filenames,
              std::string dataSet,
               bool errorIfOutOfBound = true);
  IdefixArray1D<int> dimensions;
  IdefixArray1D<int> offset;      // Actually sum_(n-1) (dimensions)
  IdefixArray1D<real> xin;

  IdefixArray1D<real> data;

  bool errorIfOutOfBound{true};

  // Fetch function that should be called inside idefix_loop
  KOKKOS_INLINE_FUNCTION
  real Get(const real x[kDim] ) const {
    int idx[kDim];
    real delta[kDim];

    for(int n = 0 ; n < kDim ; n++) {
      real xstart = xin(offset(n));
      real xend = xin(offset(n)+dimensions(n)-1);
      real x_n = x[n];

      if(std::isnan(x_n)) return(NAN);

       // Check that we're within bounds
      if(x_n < xin(offset(n))) {
        if(errorIfOutOfBound) {
          Kokkos::abort("LookupTable:: ERROR! Attempt to interpolate below your lower bound.");
        } else {
          x_n = xin(offset(n));
        }
      }

      if( x_n >= xin(offset(n)+dimensions(n)-1)) {
        if(errorIfOutOfBound) {
          Kokkos::abort("LookupTable:: ERROR! Attempt to interpolate above your upper bound.");
        } else {
          x_n = xin(offset(n)+dimensions(n)-1)*(1 - SMALL_NUMBER);
        }
      }

      // Compute index of closest element assuming even distribution
      int i = static_cast<int> ( (x_n - xstart) / (xend - xstart) * (dimensions(n)-1));

      // Check if resulting bounding elements are correct
      if(xin(offset(n) + i) > x_n || xin(offset(n) + i+1) < x_n) {
        // Nop, so the points are not evenly distributed
        // Search for the correct index (a dicotomy would be more appropriate...)
        i = 0;
        while(xin(offset(n) + i) < x_n) {
          i++;;
        }
        i = i-1; // i is overestimated by one
      }

      // Store the index
      idx[n] = i;

      // Store the elementary ratio
      delta[n] = (x_n - xin(offset(n) + i) ) / (xin(offset(n) + i+1) - xin(offset(n) + i));
    }

    // De a linear interpolation from the neightbouring points to get our value.
    real value = 0;

    // loop on all of the vertices of the neighbours
    for(unsigned int n = 0 ; n < (1 << kDim) ; n++) {
      int index = 0;
      real weight = 1.0;
      for(unsigned int m = 0 ; m < kDim ; m++) {
        index = index * dimensions(m);
        unsigned int myBit = 1 << m;
        // If bit is set, we're doing the right vertex, otherwise we're doing the left vertex
        if((n & myBit) > 0) {
          // We're on the right
          weight = weight*delta[m];
          index += idx[m]+1;
        } else {
          // We're on the left
          weight = weight*(1-delta[m]);
          index += idx[m];
        }
      }
      value = value + weight*data(index);
    }

    return(value);
  }
};

template <int kDim>
LookupTable<kDim>::LookupTable(std::vector<std::string> filenames,
                               std::string dataSet,
                               bool errOOB) {
  idfx::pushRegion("LookupTable::LookupTable");
  this->errorIfOutOfBound = errOOB;

  std::vector<uint64_t> shape;
  bool fortran_order;
  std::vector<double> dataVector;
  if(filenames.size() != kDim) {
    IDEFIX_ERROR("The list of coordinate files should match the number"
                  " of dimensions of LookupTable");
  }
  // Load the full dataset
  try {
    npy::LoadArrayFromNumpy(dataSet, shape, fortran_order, dataVector);
  } catch(std::exception &e) {
    std::stringstream errmsg;
    errmsg << e.what();
    errmsg << "LookupTable cannot load the file " << dataSet << std::endl;
    IDEFIX_ERROR(errmsg);
  }

  if(shape.size() != kDim) {
    IDEFIX_ERROR("The input numpy dataSet dimensions and LookupTable dimensions do not match");
  }
  if(fortran_order) {
    IDEFIX_ERROR("The input numpy dataSet should follow C ordering convention (not FORTRAN)");
  }

  // Load this crap in memory
  int64_t sizeTotal = 0;
  for(int n=0 ; n < shape.size() ; n++) {
    sizeTotal += shape[n];
  }

  // Allocate the required memory
  //Allocate arrays so that the data fits in it
  this->xin = IdefixArray1D<real> ("Table_x", sizeTotal);
  this->dimensions = IdefixArray1D<int> ("Table_dim", kDim);
  this->offset = IdefixArray1D<int> ("Table_offset", kDim);
  this->data =  IdefixArray1D<real> ("Table_data", dataVector.size());

  IdefixArray1D<real>::HostMirror xHost = Kokkos::create_mirror_view(this->xin);
  IdefixArray1D<int>::HostMirror dimensionsHost = Kokkos::create_mirror_view(this->dimensions);
  IdefixArray1D<int>::HostMirror offsetHost = Kokkos::create_mirror_view(this->offset);
  IdefixArray1D<real>::HostMirror dataHost = Kokkos::create_mirror_view(this->data);

  // Copy data in memory
  for(uint64_t i = 0 ; i < dataVector.size() ; i++) {
    dataHost(i) = dataVector[i];
    if(std::isnan(dataHost(i))) {
      std::stringstream msg;
      msg << "Nans were found while reading " << dataSet << std::endl;
      IDEFIX_ERROR(msg);
    }
  }

  // Copy shape arrays and coordinates
  offsetHost(0) = 0;
  for(int n = 0 ; n < kDim ; n++) {
    dimensionsHost(n) = shape[n];
    if(n>0) offsetHost(n) = offsetHost(n-1) + shape[n-1];
    std::vector<uint64_t> shapeX;
    std::vector<double> dataX;
    shapeX.clear();
    dataX.clear();
    try {
      npy::LoadArrayFromNumpy(filenames[n], shapeX, fortran_order, dataX);
    } catch(std::exception &e) {
      std::stringstream errmsg;
      errmsg << e.what() << std::endl;
      errmsg << "LookupTable cannot load the file " << filenames[n] << std::endl;
      IDEFIX_ERROR(errmsg);
    }
    if(shapeX[0] != dimensionsHost(n)) {
      idfx::cout << "ERROR: Dimension of " << filenames[n]
                 << " does not match "<< n+1 << "th dimension of " << dataSet << std::endl;
      IDEFIX_ERROR("Cannot make a lookup table out of provided numpy files");
    }
    if(fortran_order) {
      IDEFIX_ERROR("The input numpy coordinates should follow C ordering convention (not FORTRAN)");
    }
    for(int i = 0 ; i < shapeX[0] ; i++) {
      xHost(offsetHost(n)+i) = dataX[i];
      if(std::isnan(dataX[i])) {
        std::stringstream msg;
        msg << "Nans were found while reading " << filenames[n] << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
  }

  // Copy to target
  Kokkos::deep_copy(this->xin ,xHost);
  Kokkos::deep_copy(this->dimensions, dimensionsHost);
  Kokkos::deep_copy(this->offset, offsetHost);
  Kokkos::deep_copy(this->data, dataHost);

  idfx::popRegion();
}


// Constructor from CSV file
template <int kDim>
LookupTable<kDim>::LookupTable(std::string filename, char delimiter, bool errOOB) {
  idfx::pushRegion("LookupTable::LookupTable");
    this->errorIfOutOfBound = errOOB;
  if(kDim>2) {
    IDEFIX_ERROR("CSV files are only compatible with 1D and 2D tables");
  }
  // Only 1 process loads the file
  // Size of the array
  int size[2];
  // Containers for the dataset
  std::vector<real> xVector;
  std::vector<real> yVector;
  std::vector<std::vector<real>> dataVector;

  if(idfx::prank == 0) {
    std::ifstream file(filename);


    if(file.is_open()) {
      std::string line, lineWithComments;
      bool firstLine = true;
      int nx = -1;

      while(std::getline(file, lineWithComments)) {
        // get rid of comments (starting with #)
        line = lineWithComments.substr(0, lineWithComments.find("#",0));
        if (line.empty()) continue;     // skip blank line
        char firstChar = line.find_first_not_of(" ");
        if (firstChar == std::string::npos) continue;      // line is all white space
        // Walk the line
        bool firstColumn=true;
        if(kDim == 1) firstColumn = false;

        std::vector<real> dataLine;
        dataLine.clear();
        // make the line a string stream, and get all of the values separated by a delimiter
        std::stringstream str(line);
        std::string valueString;
        while(std::getline(str, valueString, delimiter)) {
          real value;
          try {
            value = std::stod(valueString);
          } catch(const std::exception& e) {
            std::stringstream errmsg;
            errmsg << e.what() << std::endl
                   << "LookupTable: Error while parsing " << filename  << ", \"" << valueString
                   << "\" cannot be converted to real." << std::endl;
            IDEFIX_ERROR(errmsg);
          }
          if(firstLine) {
            xVector.push_back(value);
          } else if(firstColumn) {
            yVector.push_back(value);
            firstColumn = false;
          } else {
            dataLine.push_back(value);
          }
        }
        // We have finished the line
        if(firstLine) {
          nx = xVector.size();
          firstLine=false;
        } else {
          if(dataLine.size() != nx) {
            IDEFIX_ERROR("LookupTable: The number of columns in the input CSV "
                          "file should be constant");
          }
          dataVector.push_back(dataLine);
          firstLine = false;
          if(kDim < 2) break; // Stop reading what's after the first two lines
        }
      }
      file.close();
      // End of file reached
    } else {
      std::stringstream errmsg;
      errmsg << "LookupTable: Unable to open file " << filename << std::endl;
      IDEFIX_ERROR(errmsg);
    }

    size[0] = xVector.size();
    if(kDim>1) {
      size[1] = yVector.size();
    } else {
      size[1] = 1;
    }
  }

  #ifdef WITH_MPI
    // Share the size of the arrays
    MPI_Bcast(size, 2, MPI_INT, 0, MPI_COMM_WORLD);
  #endif
  int sizeTotal = size[0];
  if(kDim>1) sizeTotal += size[1];

  //Allocate arrays so that the data fits in it
  this->xin = IdefixArray1D<real> ("Table_x", sizeTotal);
  this->dimensions = IdefixArray1D<int> ("Table_dim", kDim);
  this->offset = IdefixArray1D<int> ("Table_offset", kDim);
  this->data =  IdefixArray1D<real> ("Table_data", size[0]*size[1]);

  IdefixArray1D<real>::HostMirror xHost = Kokkos::create_mirror_view(this->xin);
  IdefixArray1D<int>::HostMirror dimensionsHost = Kokkos::create_mirror_view(this->dimensions);
  IdefixArray1D<int>::HostMirror offsetHost = Kokkos::create_mirror_view(this->offset);
  IdefixArray1D<real>::HostMirror dataHost = Kokkos::create_mirror_view(this->data);

  // Fill the arrays with the std::vector content
  if(idfx::prank == 0) {
    dimensionsHost(0) = size[0];
    offsetHost(0) = 0;

    for(int i = 0 ; i < xVector.size(); i++) {
      xHost(i) = xVector[i];
      if(std::isnan(xHost(i))) {
        std::stringstream msg;
        msg << "Nans were found in coordinates while reading " << filename << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
    if(kDim>1) {
      dimensionsHost(1) = size[1];
      offsetHost(1) = offsetHost(0)+dimensionsHost(0);
      for(int i = 0 ; i < yVector.size(); i++) {
        xHost(offsetHost(1)+i) = yVector[i];
        if(std::isnan(yVector[i])) {
          std::stringstream msg;
          msg << "Nans were found in coordinates while reading " << filename << std::endl;
          IDEFIX_ERROR(msg);
        }
      }
    }

    for(int j = 0 ; j < dataVector.size(); j++) {
      auto line = dataVector[j];
      for(int i = 0 ; i < line.size(); i++) {
        dataHost(i*size[1]+j) = line[i];
        if(std::isnan(line[i])) {
          std::stringstream msg;
          msg << "Nans were found in dataset while reading " << filename << std::endl;
          IDEFIX_ERROR(msg);
        }
      }
    }
  }

  #ifdef WITH_MPI
    // Share with the others
    MPI_Bcast(xHost.data(), xHost.extent(0), realMPI, 0, MPI_COMM_WORLD);
    MPI_Bcast(dimensionsHost.data(), dimensionsHost.extent(0), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offsetHost.data(), offsetHost.extent(0), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(dataHost.data(),dataHost.extent(0), realMPI, 0, MPI_COMM_WORLD);
  #endif

  // Copy to target
  Kokkos::deep_copy(this->xin ,xHost);
  Kokkos::deep_copy(this->dimensions, dimensionsHost);
  Kokkos::deep_copy(this->offset, offsetHost);
  Kokkos::deep_copy(this->data, dataHost);

  // Show the content
  /*
  idfx::cout << "x:" << std::endl;
  for(int i = 0; i < xHost.extent(0); i++) {
    idfx::cout << xHost(i) << "\t";
  }
  idfx::cout << std::endl << "y:" << std::endl;
  for(int i = 0; i < yHost.extent(0); i++) {
    idfx::cout << yHost(i) << "\t";
  }
  idfx::cout << std::endl << "data:" << std::endl;
  for(int i = 0; i < dataHost.extent(0); i++) {
    for(int j = 0; j < dataHost.extent(1); j++) {
        idfx::cout << dataHost(i,j) << "\t";
    }
    idfx::cout << std::endl;
  }*/
  idfx::popRegion();
}


#endif //UTILS_LOOKUPTABLE_HPP_
