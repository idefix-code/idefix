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

// Default transformation (and its inverse) used when the lookup table interpolates in function
// space: the interpolation is then performed on log(x) and log(data), and the interpolated value
// is transformed back with exp.
// Any functor exposing a KOKKOS_INLINE_FUNCTION operator()(const real) can be used instead, so
// that the transformation is available both on the host and on the device.
struct LookupTableLog {
  KOKKOS_INLINE_FUNCTION real operator() (const real x) const {
    return(LOG(x));
  }
};

struct LookupTableExp {
  KOKKOS_INLINE_FUNCTION real operator() (const real x) const {
    return(EXP(x));
  }
};

// Neighbours of a point of a lookup table, as they are found by the search performed by Get:
// idx[n] is the index of the left neighbour along the dimension n, and delta[n] is the elementary
// ratio used to weight the two neighbours of that dimension.
// Get (and GetHost) fill this structure when it is given as an argument, and GetNeighbours and
// GetNeighboursIndx then reuse the search it contains instead of performing it a second time
// (which is only done when they are called with the same coordinates: they search the table again
// as usual when a different x is requested).
// This structure is kept by the caller, so that it is thread-private when it is declared inside
// an idefix_for loop. The lookup table itself is shared by all of the threads of a loop, and can
// therefore not be used to store anything of the sort.
template <const int kDim>
struct LookupTableNeighbours {
  real x[kDim];        // coordinates for which the neighbours below were computed
  int idx[kDim];       // index of the left neighbour along each dimension
  real delta[kDim];    // elementary ratio between the two neighbours of each dimension
  bool valid{false};   // whether the neighbours above have been successfully computed

  // Check whether this structure already holds the neighbours of the coordinates xIn
  KOKKOS_INLINE_FUNCTION
  bool Matches(const real xIn[kDim]) const {
    if(!valid) return(false);
    for(int n = 0 ; n < kDim ; n++) {
      if(x[n] != xIn[n]) return(false);
    }
    return(true);
  }
};

template <const int kDim, class TFunc = LookupTableLog, class TInvFunc = LookupTableExp>
class LookupTable {
 public:
  LookupTable() = default;
  // Constructor from a CSV/ASCII file.
  // By default (readColumns=false), the arrays are read as *lines* of the input file:
  //    1D: 1st line = coordinates (xinHost), 2nd line = data (dataHost)
  //    2D: 1st line = coordinates of the 1st dimension, then each line starts with the
  //        coordinate of the 2nd dimension, followed by the data of that line.
  // For 1D tables only, the arrays can also be read as *columns* of the input file by setting
  // readColumns=true:
  //    1st column = coordinates (xinHost), 2nd column = data (dataHost)
  // All of the constructors accept the optional argument interpolateInFuncSpace: when it is set
  // to true, the coordinates and the data are stored as func(x) and func(data), the interpolation
  // is performed in that space, and Get returns invFunc of the interpolated value (see below).
  LookupTable(std::string filename, char delimiter, bool errorIfOutOfBound = true,
              bool readColumns = false, bool interpolateInFuncSpace = false,
              TFunc func = TFunc(), TInvFunc invFunc = TInvFunc());
  LookupTable(std::vector<std::string> filenames,
              std::string dataSet,
               bool errorIfOutOfBound = true, bool interpolateInFuncSpace = false,
               TFunc func = TFunc(), TInvFunc invFunc = TInvFunc());
  template<typename T, typename ... Args>
  LookupTable(Kokkos::View<T, Args...> array,
              std::array<IdefixHostArray1D<real>,kDim>,
               bool errorIfOutOfBound = true, bool interpolateInFuncSpace = false,
               TFunc func = TFunc(), TInvFunc invFunc = TInvFunc());

  IdefixArray1D<int> dimensionsDev;
  IdefixArray1D<int> offsetDev;      // Actually sum_(n-1) (dimensions)
  IdefixArray1D<real> xinDev;
  IdefixArray1D<real> dataDev;

  IdefixHostArray1D<int> dimensionsHost;
  IdefixHostArray1D<int> offsetHost;      // Actually sum_(n-1) (dimensions)
  IdefixHostArray1D<real> xinHost;
  IdefixHostArray1D<real> dataHost;


  bool errorIfOutOfBound{true};

  // When enabled, the table is interpolated in function space: xin and data store func(x) and
  // func(data), and Get returns invFunc(interpolated value). func and invFunc are stored by
  // value, so that they are available on the device as well as on the host.
  bool interpolateInFuncSpace{false};
  TFunc func;
  TInvFunc invFunc;

  // Transform the coordinates and the data of the table in function space.
  // This is called on the host by the constructors, before the arrays are copied on the device.
  void ToFuncSpace() {
    for(int i = 0 ; i < xinHost.extent(0) ; i++) {
      xinHost(i) = func(xinHost(i));
      if(!std::isfinite(xinHost(i))) {
        IDEFIX_ERROR("LookupTable: the transformation of the coordinates in function space "
                     "produced invalid values (with the default log, the coordinates of the "
                     "table should all be strictly positive)");
      }
    }
    for(int i = 0 ; i < dataHost.extent(0) ; i++) {
      dataHost(i) = func(dataHost(i));
      if(!std::isfinite(dataHost(i))) {
        IDEFIX_ERROR("LookupTable: the transformation of the data in function space "
                     "produced invalid values (with the default log, the data of the "
                     "table should all be strictly positive)");
      }
    }
  }

  // Generic search of the neighbours used for the interpolation, for all kinds of input arrays.
  // On output, idx[n] is the index of the left neighbour along the dimension n (so that the
  // interpolation is performed between xin(offset(n)+idx[n]) and xin(offset(n)+idx[n]+1)), and
  // delta[n] is the elementary ratio used to weight these two neighbours.
  // Returns false when the requested coordinates contain nans, and true otherwise.
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  bool GetIndices(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin,
                  int idx[kDim], real delta[kDim]) const {
    for(int n = 0 ; n < kDim ; n++) {
      real xstart = xin(offset(n));
      real xend = xin(offset(n)+dimensions(n)-1);
      // When interpolating in function space, xin already stores func(x), so the coordinate
      // we are looking for should be transformed accordingly
      real x_n = interpolateInFuncSpace ? func(x[n]) : x[n];

      if(std::isnan(x_n)) return(false);

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

    return(true);
  }

  // Index in the data array of the vertex "vertex" of the neighbours. Each bit of "vertex" tells
  // whether we are on the right (1) or on the left (0) of the corresponding dimension, so that
  // there are 2^kDim vertices surrounding the point we are interpolating.
  template<typename Tint>
  KOKKOS_INLINE_FUNCTION
  int GetDataIndex(Tint &dimensions, const int idx[kDim], const unsigned int vertex) const {
    int index = 0;
    for(unsigned int m = 0 ; m < kDim ; m++) {
      index = index * dimensions(m);
      unsigned int myBit = 1 << m;
      // If bit is set, we're doing the right vertex, otherwise we're doing the left vertex
      if((vertex & myBit) > 0) {
        index += idx[m]+1;
      } else {
        index += idx[m];
      }
    }
    return(index);
  }

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
      // When interpolating in function space, xin already stores func(x), so the coordinate
      // we are looking for should be transformed accordingly
      real x_n = interpolateInFuncSpace ? func(x[n]) : x[n];

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

    // Do a linear interpolation from the neightbouring points to get our value.
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

    // The interpolation was performed on func(data), so we transform the result back
    if(interpolateInFuncSpace) value = invFunc(value);

    return(value);
  }

  // Fill "neighbours" with the neighbours of x, unless it already holds them (in which case the
  // table is not searched again). This is what allows Get, GetNeighbours and GetNeighboursIndx to
  // share a single search when they are called successively with the same coordinates.
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  void SearchNeighbours(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin,
                        LookupTableNeighbours<kDim> &neighbours) const {
    // Nothing to do if the neighbours of these very coordinates are already known
    if(neighbours.Matches(x)) return;

    neighbours.valid = GetIndices(x, dimensions, offset, xin, neighbours.idx, neighbours.delta);
    for(int n = 0 ; n < kDim ; n++) {
      neighbours.x[n] = x[n];
    }
  }

  // Generic getter which stores in "neighbours" the elements of the table it used, so that a
  // subsequent call to GetNeighbours or GetNeighboursIndx with the same coordinates does not
  // search the table again
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  real Get(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin, Treal &data,
           LookupTableNeighbours<kDim> &neighbours) const {
    SearchNeighbours(x, dimensions, offset, xin, neighbours);

    if(!neighbours.valid) return(NAN);

    // Do a linear interpolation from the neightbouring points to get our value.
    real value = 0;

    // loop on all of the vertices of the neighbours
    for(unsigned int n = 0 ; n < (1 << kDim) ; n++) {
      real weight = 1.0;
      for(unsigned int m = 0 ; m < kDim ; m++) {
        unsigned int myBit = 1 << m;
        // If bit is set, we're doing the right vertex, otherwise we're doing the left vertex
        if((n & myBit) > 0) {
          // We're on the right
          weight = weight*neighbours.delta[m];
        } else {
          // We're on the left
          weight = weight*(1-neighbours.delta[m]);
        }
      }
      value = value + weight*data(GetDataIndex(dimensions, neighbours.idx, n));
    }

    // The interpolation was performed on func(data), so we transform the result back
    if(interpolateInFuncSpace) value = invFunc(value);

    return(value);
  }

  // Generic getter for the neighbours used by the interpolation, for all kinds of input arrays.
  // On output, xN[2*n] and xN[2*n+1] are the coordinates bracketing x[n] along the dimension n,
  // and dataN[v] is the data at the vertex v of these neighbours (see GetDataIndex for the
  // convention on v). When the table is interpolated in function space, both are returned in the
  // original space of the table (i.e. invFunc is applied to the values stored in the table).
  // All of the outputs are set to nan when the requested coordinates contain nans.
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  void GetNeighbours(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin, Treal &data,
                     real xN[2*kDim], real dataN[1 << kDim]) const {
    LookupTableNeighbours<kDim> neighbours;
    GetNeighbours(x, dimensions, offset, xin, data, neighbours, xN, dataN);
  }

  // Same as above, but the search is stored in (and reused from) "neighbours": the table is only
  // searched again when "neighbours" does not already hold the neighbours of x, e.g. because it
  // was filled by a previous call to Get with these very same coordinates.
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  void GetNeighbours(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin, Treal &data,
                     LookupTableNeighbours<kDim> &neighbours,
                     real xN[2*kDim], real dataN[1 << kDim]) const {
    SearchNeighbours(x, dimensions, offset, xin, neighbours);

    if(!neighbours.valid) {
      for(int n = 0 ; n < 2*kDim ; n++) xN[n] = NAN;
      for(unsigned int n = 0 ; n < (1 << kDim) ; n++) dataN[n] = NAN;
      return;
    }

    // Coordinates of the neighbours along each dimension
    for(int n = 0 ; n < kDim ; n++) {
      xN[2*n] = xin(offset(n) + neighbours.idx[n]);
      xN[2*n+1] = xin(offset(n) + neighbours.idx[n]+1);
      if(interpolateInFuncSpace) {
        xN[2*n] = invFunc(xN[2*n]);
        xN[2*n+1] = invFunc(xN[2*n+1]);
      }
    }

    // Data on each vertex of the neighbours
    for(unsigned int n = 0 ; n < (1 << kDim) ; n++) {
      dataN[n] = data(GetDataIndex(dimensions, neighbours.idx, n));
      if(interpolateInFuncSpace) dataN[n] = invFunc(dataN[n]);
    }
  }

  // Generic getter for the indices of the neighbours used by the interpolation, for all kinds of
  // input arrays. On output, idx[n] is the index of the left neighbour along the dimension n
  // (the right one being idx[n]+1), and dataIdx[v] is the index in the data array of the vertex v
  // of these neighbours (see GetDataIndex for the convention on v).
  // All of the outputs are set to -1 when the requested coordinates contain nans.
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndx(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin,
                         int idx[kDim], int dataIdx[1 << kDim]) const {
    LookupTableNeighbours<kDim> neighbours;
    GetNeighboursIndx(x, dimensions, offset, xin, neighbours, idx, dataIdx);
  }

  // Same as above, but the search is stored in (and reused from) "neighbours", exactly like the
  // GetNeighbours variant above
  template<typename Tint, typename Treal>
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndx(const real x[kDim], Tint &dimensions, Tint &offset, Treal &xin,
                         LookupTableNeighbours<kDim> &neighbours,
                         int idx[kDim], int dataIdx[1 << kDim]) const {
    SearchNeighbours(x, dimensions, offset, xin, neighbours);

    if(!neighbours.valid) {
      for(int n = 0 ; n < kDim ; n++) idx[n] = -1;
      for(unsigned int n = 0 ; n < (1 << kDim) ; n++) dataIdx[n] = -1;
      return;
    }

    for(int n = 0 ; n < kDim ; n++) {
      idx[n] = neighbours.idx[n];
    }
    for(unsigned int n = 0 ; n < (1 << kDim) ; n++) {
      dataIdx[n] = GetDataIndex(dimensions, neighbours.idx, n);
    }
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

  // Getter for the neighbours used by the interpolation, on device
  KOKKOS_INLINE_FUNCTION
  void GetNeighbours(const real x[kDim], real xN[2*kDim], real dataN[1 << kDim]) const {
    GetNeighbours(x, dimensionsDev, offsetDev, xinDev, dataDev, xN, dataN);
  }

  // Getter for the neighbours used by the interpolation, on Host
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursHost(const real x[kDim], real xN[2*kDim], real dataN[1 << kDim]) const {
    GetNeighbours(x, dimensionsHost, offsetHost, xinHost, dataHost, xN, dataN);
  }

  // Getter for the indices of the neighbours used by the interpolation, on device
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndx(const real x[kDim], int idx[kDim], int dataIdx[1 << kDim]) const {
    GetNeighboursIndx(x, dimensionsDev, offsetDev, xinDev, idx, dataIdx);
  }

  // Getter for the indices of the neighbours used by the interpolation, on Host
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndxHost(const real x[kDim], int idx[kDim], int dataIdx[1 << kDim]) const {
    GetNeighboursIndx(x, dimensionsHost, offsetHost, xinHost, idx, dataIdx);
  }

  // Getter on device, which stores the neighbours it used in "neighbours". Giving that same
  // structure to GetNeighbours or GetNeighboursIndx below then avoids searching the table twice.
  KOKKOS_INLINE_FUNCTION
  real Get(const real x[kDim], LookupTableNeighbours<kDim> &neighbours) const {
    return(Get(x, dimensionsDev, offsetDev, xinDev, dataDev, neighbours));
  }

  // Getter on Host, which stores the neighbours it used in "neighbours"
  KOKKOS_INLINE_FUNCTION
  real GetHost(const real x[kDim], LookupTableNeighbours<kDim> &neighbours) const {
    return(Get(x, dimensionsHost, offsetHost, xinHost, dataHost, neighbours));
  }

  // Getter for the neighbours used by the interpolation, on device, reusing the search stored in
  // "neighbours" when it was performed for these very same coordinates
  KOKKOS_INLINE_FUNCTION
  void GetNeighbours(const real x[kDim], LookupTableNeighbours<kDim> &neighbours,
                     real xN[2*kDim], real dataN[1 << kDim]) const {
    GetNeighbours(x, dimensionsDev, offsetDev, xinDev, dataDev, neighbours, xN, dataN);
  }

  // Getter for the neighbours used by the interpolation, on Host, reusing the search stored in
  // "neighbours" when it was performed for these very same coordinates
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursHost(const real x[kDim], LookupTableNeighbours<kDim> &neighbours,
                         real xN[2*kDim], real dataN[1 << kDim]) const {
    GetNeighbours(x, dimensionsHost, offsetHost, xinHost, dataHost, neighbours, xN, dataN);
  }

  // Getter for the indices of the neighbours, on device, reusing the search stored in
  // "neighbours" when it was performed for these very same coordinates
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndx(const real x[kDim], LookupTableNeighbours<kDim> &neighbours,
                         int idx[kDim], int dataIdx[1 << kDim]) const {
    GetNeighboursIndx(x, dimensionsDev, offsetDev, xinDev, neighbours, idx, dataIdx);
  }

  // Getter for the indices of the neighbours, on Host, reusing the search stored in
  // "neighbours" when it was performed for these very same coordinates
  KOKKOS_INLINE_FUNCTION
  void GetNeighboursIndxHost(const real x[kDim], LookupTableNeighbours<kDim> &neighbours,
                             int idx[kDim], int dataIdx[1 << kDim]) const {
    GetNeighboursIndx(x, dimensionsHost, offsetHost, xinHost, neighbours, idx, dataIdx);
  }
};

template <int kDim, class TFunc, class TInvFunc>
LookupTable<kDim, TFunc, TInvFunc>::LookupTable(std::vector<std::string> filenames,
                               std::string dataSet,
                               bool errOOB, bool funcSpace, TFunc funcIn, TInvFunc invFuncIn) {
  idfx::pushRegion("LookupTable::LookupTable");
  this->errorIfOutOfBound = errOOB;
  this->interpolateInFuncSpace = funcSpace;
  this->func = funcIn;
  this->invFunc = invFuncIn;

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
  // Allocate arrays so that the data fits in it
  this->xinDev = IdefixArray1D<real> ("Table_x", sizeTotal);
  this->dimensionsDev = IdefixArray1D<int> ("Table_dim", kDim);
  this->offsetDev = IdefixArray1D<int> ("Table_offset", kDim);
  this->dataDev =  IdefixArray1D<real> ("Table_data", dataVector.size());

  this->xinHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dataDev);

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

  // Transform the table in function space if required
  if(this->interpolateInFuncSpace) this->ToFuncSpace();

  // Copy to target
  Kokkos::deep_copy(this->xinDev ,xinHost);
  Kokkos::deep_copy(this->dimensionsDev, dimensionsHost);
  Kokkos::deep_copy(this->offsetDev, offsetHost);
  Kokkos::deep_copy(this->dataDev, dataHost);

  idfx::popRegion();
}


// Constructor from CSV file
// The coordinates and the data are read as lines of the input file, unless readColumns is set
// to true, in which case they are read as columns of the input file (1D tables only, see the
// declaration of the class for a description of both layouts).
template <int kDim, class TFunc, class TInvFunc>
LookupTable<kDim, TFunc, TInvFunc>::LookupTable(std::string filename, char delimiter, bool errOOB,
                               bool readColumns, bool funcSpace, TFunc funcIn, TInvFunc invFuncIn) {
  idfx::pushRegion("LookupTable::LookupTable");
    this->errorIfOutOfBound = errOOB;
    this->interpolateInFuncSpace = funcSpace;
    this->func = funcIn;
    this->invFunc = invFuncIn;
  if(kDim>2) {
    IDEFIX_ERROR("CSV files are only compatible with 1D and 2D tables");
  }
  if(kDim>1 && readColumns) {
    IDEFIX_ERROR("CSV files can only be read as columns for 1D tables");
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
      // Full content of the file, stored as a list of lines, each line being a list of values
      std::vector<std::vector<real>> fileContent;

      while(std::getline(file, lineWithComments)) {
        // get rid of comments (starting with #)
        line = lineWithComments.substr(0, lineWithComments.find("#",0));
        if (line.empty()) continue;     // skip blank line
        std::size_t firstChar = line.find_first_not_of(" ");
        if (firstChar == std::string::npos) continue;      // line is all white space
        // Walk the line
        std::vector<real> lineVector;
        lineVector.clear();
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
          lineVector.push_back(value);
        }
        // We have finished the line
        fileContent.push_back(lineVector);
        // When read as lines, a 1D table is fully described by the first two lines of the file,
        // so we stop reading what's after them
        if(kDim < 2 && !readColumns && fileContent.size() == 2) break;
      }
      file.close();
      // End of file reached

      if(fileContent.size() < 2) {
        IDEFIX_ERROR("LookupTable: The input CSV file should contain at least two lines");
      }

      // Dispatch the content of the file in the coordinate and data containers.
      // Note that dataVector is always indexed as dataVector[j][i], where i (resp. j) is the
      // index along the 1st (resp. 2nd) dimension of the table.
      if(readColumns) {
        // (1D tables only) the coordinates are stored in the 1st column of the input file,
        // while the data is stored in its 2nd column
        dataVector.push_back(std::vector<real>());
        for(int i = 0 ; i < fileContent.size() ; i++) {
          if(fileContent[i].size() < 2) {
            IDEFIX_ERROR("LookupTable: The input CSV file should have at least two columns "
                          "when the table is read as columns");
          }
          xVector.push_back(fileContent[i][0]);
          dataVector[0].push_back(fileContent[i][1]);
        }
      } else {
        // (default) the arrays are stored as lines of the input file
        // 1st line always gives the coordinates of the 1st dimension
        xVector = fileContent[0];
        const int nx = xVector.size();
        if(kDim < 2) {
          // 2nd line is the data
          if(fileContent[1].size() != nx) {
            IDEFIX_ERROR("LookupTable: The number of columns in the input CSV "
                          "file should be constant");
          }
          dataVector.push_back(fileContent[1]);
        } else {
          // each of the following lines starts with the coordinate of the 2nd dimension,
          // followed by the data of that line
          for(int j = 1 ; j < fileContent.size() ; j++) {
            if(fileContent[j].size() != nx+1) {
              IDEFIX_ERROR("LookupTable: The number of columns in the input CSV "
                            "file should be constant");
            }
            yVector.push_back(fileContent[j][0]);
            dataVector.push_back(std::vector<real>(fileContent[j].begin()+1,
                                                   fileContent[j].end()));
          }
        }
      }
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

  this->xinHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dataDev);

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

  // Transform the table in function space if required
  if(this->interpolateInFuncSpace) this->ToFuncSpace();

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
template<const int kDim, class TFunc, class TInvFunc>
template<typename T, typename ... Args>
LookupTable<kDim, TFunc, TInvFunc>::LookupTable(Kokkos::View<T, Args...> array,
            std::array<IdefixHostArray1D<real>,kDim> x,
              bool errOOB, bool funcSpace, TFunc funcIn, TInvFunc invFuncIn) {
  idfx::pushRegion("LookupTable::LookupTable");
  this->errorIfOutOfBound = errOOB;
  this->interpolateInFuncSpace = funcSpace;
  this->func = funcIn;
  this->invFunc = invFuncIn;

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

  this->xinHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->xinDev);
  this->dimensionsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dimensionsDev);
  this->offsetHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->offsetDev);
  this->dataHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), this->dataDev);

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

  // Transform the table in function space if required
  if(this->interpolateInFuncSpace) this->ToFuncSpace();

  // Copy to target
  Kokkos::deep_copy(this->xinDev ,xinHost);
  Kokkos::deep_copy(this->dimensionsDev, dimensionsHost);
  Kokkos::deep_copy(this->offsetDev, offsetHost);
  Kokkos::deep_copy(this->dataDev, dataHost);

  idfx::popRegion();
              }

#endif //UTILS_LOOKUPTABLE_HPP_
