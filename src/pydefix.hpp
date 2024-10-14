// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef PYDEFIX_HPP_
#define PYDEFIX_HPP_


#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include "idefix.hpp"
#include "input.hpp"

namespace py = pybind11;

class DataBlock;
class DataBlockHost;

class Pydefix {
 public:
  explicit Pydefix(Input&);
  ~Pydefix();
  void Output(DataBlock &, int);
  void InitFlow(DataBlock &);
  void ShowConfig();
  bool isActive{false};
  bool haveOutput{false};
  bool haveInitflow{false};
 private:
  template<typename... Ts>
  void CallScript(std::string, std::string, Ts...);
  static int ninstance;
  std::string scriptFilename;
  std::string outputFunctionName;
  std::string initflowFunctionName;
};


namespace pybind11 { namespace detail {
// Caster for IdefixArray4D<T>
template <typename T> struct type_caster<IdefixHostArray4D<T>> {
 public:
  PYBIND11_TYPE_CASTER(IdefixHostArray4D<T>, _("IdefixHostArray4D<T>"));

  // Conversion part 1 (Python -> C++)
  bool load(py::handle src, bool convert) {
    if ( !convert && !py::array_t<T>::check_(src) )
      return false;

    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    if ( !buf )
      return false;

    auto dims = buf.ndim();
    if ( dims != 4  )
      return false;

    std::vector<size_t> shape(4);

    for ( int i = 0 ; i < 4 ; ++i )
      shape[i] = buf.shape()[i];


    value = IdefixHostArray4D<T>("pyArray",shape[0], shape[1], shape[2], shape[3]);

    // Still need to fill in with buf.data()+buf.size()
    IDEFIX_ERROR("Python->Idefix Not implemented");

    return true;
  }

  //Conversion part 2 (C++ -> Python)
  static py::handle cast(const IdefixHostArray4D<T>& src,
                         py::return_value_policy policy,
                         py::handle parent) {
    py::none dummyDataOwner;
    py::array_t<real, py::array::c_style> a({src.extent(0),
                                             src.extent(1),
                                             src.extent(2),
                                             src.extent(3)},
                                             src.data(), dummyDataOwner);

    return a.release();
  }
};
// Caster for IdefixArray3D<T>
template <typename T> struct type_caster<IdefixHostArray3D<T>> {
 public:
  PYBIND11_TYPE_CASTER(IdefixHostArray3D<T>, _("IdefixHostArray3D<T>"));

  // Conversion part 1 (Python -> C++)
  bool load(py::handle src, bool convert) {
    if ( !convert && !py::array_t<T>::check_(src) )
      return false;

    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    if ( !buf )
      return false;

    auto dims = buf.ndim();
    if ( dims != 3  )
      return false;

    std::vector<size_t> shape(3);

    for ( int i = 0 ; i < 3 ; ++i )
      shape[i] = buf.shape()[i];


    value = IdefixHostArray3D<T>("pyArray",shape[0], shape[1], shape[2]);

    // Still need to fill in with buf.data()+buf.size()
    IDEFIX_ERROR("Python->Idefix Not implemented");
    return true;
  }

  //Conversion part 2 (C++ -> Python)
  static py::handle cast(const IdefixHostArray3D<T>& src,
                         py::return_value_policy policy,
                         py::handle parent) {
    py::none dummyDataOwner;
    py::array_t<real, py::array::c_style> a({src.extent(0),
                                             src.extent(1),
                                             src.extent(2)},
                                             src.data(),dummyDataOwner);
    return a.release();
  }
};
// Caster for IdefixArray2D<T>
template <typename T> struct type_caster<IdefixHostArray2D<T>> {
 public:
  PYBIND11_TYPE_CASTER(IdefixHostArray2D<T>, _("IdefixHostArray2D<T>"));

  // Conversion part 1 (Python -> C++)
  bool load(py::handle src, bool convert) {
    if ( !convert && !py::array_t<T>::check_(src) )
      return false;

    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    if ( !buf )
      return false;

    auto dims = buf.ndim();
    if ( dims != 2  )
      return false;

    std::vector<size_t> shape(2);

    for ( int i = 0 ; i < 2 ; ++i )
      shape[i] = buf.shape()[i];


    value = IdefixHostArray2D<T>("pyArray",shape[0], shape[1]);

    // Still need to fill in with buf.data()+buf.size()
    IDEFIX_ERROR("Python->Idefix Not implemented");
    return true;
  }

  //Conversion part 2 (C++ -> Python)
  static py::handle cast(const IdefixHostArray2D<T>& src,
                         py::return_value_policy policy,
                         py::handle parent) {
    py::none dummyOwner;
    py::array_t<real, py::array::c_style> a({src.extent(0),src.extent(1)},src.data(),dummyOwner);
    return a.release();
  }
};
// Caster for IdefixArray1D<T>
template <typename T> struct type_caster<IdefixHostArray1D<T>> {
 public:
  PYBIND11_TYPE_CASTER(IdefixHostArray1D<T>, _("IdefixHostArray1D<T>"));

  // Conversion part 1 (Python -> C++)
  bool load(py::handle src, bool convert) {
    if ( !convert && !py::array_t<T>::check_(src) )
      return false;

    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    if ( !buf )
      return false;

    auto dims = buf.ndim();
    if ( dims != 1  )
      return false;

    std::vector<size_t> shape(1);

    for ( int i = 0 ; i < 1 ; ++i )
      shape[i] = buf.shape()[i];


    value = IdefixHostArray1D<T>("pyArray",shape[0]);

    // Still need to fill in with buf.data()+buf.size()
    IDEFIX_ERROR("Python->Idefix Not implemented");
    return true;
  }

  //Conversion part 2 (C++ -> Python)
  static py::handle cast(const IdefixHostArray1D<T>& src,
                         py::return_value_policy policy,
                         py::handle parent) {
    py::none dummyDataOwner;
    py::array_t<real, py::array::c_style> a(src.extent(0),src.data(),dummyDataOwner);
    return a.release();
  }
};
} // namespace detail
} // namespace pybind11

#endif // PYDEFIX_HPP_
