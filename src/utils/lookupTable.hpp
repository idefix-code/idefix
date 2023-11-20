// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
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
  template<typename T, typename ... Args>
  LookupTable(Kokkos::View<T, Args...> array,
              std::array<IdefixHostArray1D<real>,kDim>,
               bool errorIfOutOfBound = true);

  IdefixArray1D<int> dimensionsDev;
  IdefixArray1D<int> offsetDev;      // Actually sum_(n-1) (dimensions)
  IdefixArray1D<real> xinDev;
  IdefixArray1D<real> dataDev;

  IdefixHostArray1D<int> dimensionsHost;
  IdefixHostArray1D<int> offsetHost;      // Actually sum_(n-1) (dimensions)
  IdefixHostArray1D<real> xinHost;
  IdefixHostArray1D<real> dataHost;


  bool errorIfOutOfBound{true};

  // Generic getter for all kinds of input arrays
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  real Get(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin, Treal &data) const {
  // Fetch function that should be called inside idefix_loop
    int idx[kDim];
    real delta[kDim];

    for(int n = 0 ; n < kDim ; n++) {
      real xstart = xin(offset(n));
      real xend = xin(offset(n)+dimensions(n)-1);
      real x_n = x[n];

      if(std::isnan(x_n)) return(NAN);

      // Compute index of closest element assuming even distribution
      int i;

       // Check that we're within bounds
      if(x_n < xstart) {
        if(errorIfOutOfBound) {
          Kokkos::abort("LookupTable:: ERROR! Attempt to interpolate below your lower bound.");
        } else {
          x_n = xstart;
          i = 0;
        }
      } else if( x_n > xend) {
        if(errorIfOutOfBound) {
          Kokkos::abort("LookupTable:: ERROR! Attempt to interpolate above your upper bound.");
        } else {
          // We set x_n=xend, and we do the interpolation between xin(dim-2) and xin(dim-1),
          // so i= dim-2
          i = dimensions(n)-2;
          x_n = xend;
        }
      } else {
        // Bounds are fine,
        i = static_cast<int> ( (x_n - xstart) / (xend - xstart) * (dimensions(n)-1));
        if(i == dimensions(n)-1) i = dimensions(n)-2;
        // Check if resulting bounding elements are correct
        if(xin(offset(n) + i) > x_n || xin(offset(n) + i+1) < x_n) {
          // Nop, so the points are not evenly distributed
          // Search for the correct index (a dicotomy would be more appropriate...)

          i = 0;
          while(xin(offset(n) + i) < x_n && i < dimensions(n)-1 ) {
            i++;
          }
          i = i-1; // i is overestimated by one
        }
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

  // Getter on device
  KOKKOS_INLINE_FUNCTION
  real Get(const real x[kDim]) const {
    return(Get(x, dimensionsDev, offsetDev, xinDev, dataDev));
  }

  // Getter on Host
  KOKKOS_INLINE_FUNCTION
  real GetHost(const real x[kDim]) const {
    return(Get(x, dimensionsHost, offsetHost, xinHost, dataHost));
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
  this->xinDev = IdefixArray1D<real> ("Table_x", sizeTotal);
  this->dimensionsDev = IdefixArray1D<int> ("Table_dim", kDim);
  this->offsetDev = IdefixArray1D<int> ("Table_offset", kDim);
  this->dataDev =  IdefixArray1D<real> ("Table_data", dataVector.size());

  this->xinHost = Kokkos::create_mirror_view(this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(this->dataDev);

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
      xinHost(offsetHost(n)+i) = dataX[i];
      if(std::isnan(dataX[i])) {
        std::stringstream msg;
        msg << "Nans were found while reading " << filenames[n] << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
  }

  // Copy to target
  Kokkos::deep_copy(this->xinDev ,xinHost);
  Kokkos::deep_copy(this->dimensionsDev, dimensionsHost);
  Kokkos::deep_copy(this->offsetDev, offsetHost);
  Kokkos::deep_copy(this->dataDev, dataHost);

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
  this->xinDev = IdefixArray1D<real> ("Table_x", sizeTotal);
  this->dimensionsDev = IdefixArray1D<int> ("Table_dim", kDim);
  this->offsetDev = IdefixArray1D<int> ("Table_offset", kDim);
  this->dataDev =  IdefixArray1D<real> ("Table_data", size[0]*size[1]);

  this->xinHost = Kokkos::create_mirror_view(this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(this->dataDev);

  // Fill the arrays with the std::vector content
  if(idfx::prank == 0) {
    dimensionsHost(0) = size[0];
    offsetHost(0) = 0;

    for(int i = 0 ; i < xVector.size(); i++) {
      xinHost(i) = xVector[i];
      if(std::isnan(xinHost(i))) {
        std::stringstream msg;
        msg << "Nans were found in coordinates while reading " << filename << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
    if(kDim>1) {
      dimensionsHost(1) = size[1];
      offsetHost(1) = offsetHost(0)+dimensionsHost(0);
      for(int i = 0 ; i < yVector.size(); i++) {
        xinHost(offsetHost(1)+i) = yVector[i];
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
    MPI_Bcast(xinHost.data(), xinHost.extent(0), realMPI, 0, MPI_COMM_WORLD);
    MPI_Bcast(dimensionsHost.data(), dimensionsHost.extent(0), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offsetHost.data(), offsetHost.extent(0), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(dataHost.data(),dataHost.extent(0), realMPI, 0, MPI_COMM_WORLD);
  #endif

  // Copy to target
  Kokkos::deep_copy(this->xinDev ,xinHost);
  Kokkos::deep_copy(this->dimensionsDev, dimensionsHost);
  Kokkos::deep_copy(this->offsetDev, offsetHost);
  Kokkos::deep_copy(this->dataDev, dataHost);

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


// Constructor from IdefixHostArray
template<const int kDim>
template<typename T, typename ... Args>
LookupTable<kDim>::LookupTable(Kokkos::View<T, Args...> array,
            std::array<IdefixHostArray1D<real>,kDim> x,
              bool errOOB) {
  idfx::pushRegion("LookupTable::LookupTable");
  this->errorIfOutOfBound = errOOB;

  std::vector<uint64_t> shape(kDim);
  for(int i = 0 ; i < kDim ; i++) shape[i] = x[i].extent(0);

  // Load this crap in memory
  int64_t sizeX = 0;
  int64_t sizeTotal = 1;
  for(int n=0 ; n < shape.size() ; n++) {
    sizeX += shape[n];
    sizeTotal *= shape[n];
    if(array.extent(kDim-n-1) != shape[n]) {
      std::stringstream errmsg;
      errmsg << "The " << n+1 << "th dimension of your input array (" << array.extent(n)
             << ") does not match the size of the corresponding x vector (" << shape[n];
      IDEFIX_ERROR(errmsg);
    }
  }

  // Allocate the required memory
  //Allocate arrays so that the data fits in it
  this->xinDev = IdefixArray1D<real> ("Table_x", sizeX);
  this->dimensionsDev = IdefixArray1D<int> ("Table_dim", kDim);
  this->offsetDev = IdefixArray1D<int> ("Table_offset", kDim);
  this->dataDev =  IdefixArray1D<real> ("Table_data", sizeTotal);

  this->xinHost = Kokkos::create_mirror_view(this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(this->dataDev);

  // Copy data in memory
  for(uint64_t n = 0 ; n < sizeTotal ; n++) {
    real q;
    if constexpr(kDim == 1) {
      q = array(n);
    } else if constexpr(kDim == 2) {
      int i = n  / shape[1];
      int j = (n - i * shape[1]);
      q = array(j,i);
    } else if constexpr(kDim == 3) {
      int i = n / (shape[2]*shape[1]);
      int j = (n - i * shape[2]*shape[1]) / shape[2];
      int k = (n - i * shape[2]*shape[1] - j * shape[2]);
      q = array(k,j,i);
    } else {
      IDEFIX_ERROR("The lookup table only handles array of rank <= 3");
    }
    if(std::isnan(q)) {
      std::stringstream msg;
      msg << "Nans were found while loading the array." <<  std::endl;
      IDEFIX_ERROR(msg);
    }
    dataHost(n) = q;
  }

  // Copy shape arrays and coordinates
  offsetHost(0) = 0;
  for(int n = 0 ; n < kDim ; n++) {
    dimensionsHost(n) = shape[n];
    if(n>0) offsetHost(n) = offsetHost(n-1) + shape[n-1];

    for(int i = 0 ; i < shape[n] ; i++) {
      xinHost(offsetHost(n)+i) = x[n](i);
      if(std::isnan(x[n](i))) {
        std::stringstream msg;
        msg << "Nans were found while reading x[" << n << "]." << std::endl;
        IDEFIX_ERROR(msg);
      }
    }
  }

  // Copy to target
  Kokkos::deep_copy(this->xinDev ,xinHost);
  Kokkos::deep_copy(this->dimensionsDev, dimensionsHost);
  Kokkos::deep_copy(this->offsetDev, offsetHost);
  Kokkos::deep_copy(this->dataDev, dataHost);

  idfx::popRegion();
              }

#endif //UTILS_LOOKUPTABLE_HPP_
