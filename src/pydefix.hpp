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
#include "idefix.hpp"
#include "pydefix.hpp"

namespace py = pybind11;

class Pydefix {
 public:
  Pydefix();
  ~Pydefix();
  py::array_t<real> toNumpyArray(const IdefixHostArray3D<real>&);
  py::array_t<real> toNumpyArray(const IdefixHostArray4D<real>&);

 private:
  static int ninstance;
};

// Caster for IdefixArray4D<T>
namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<IdefixHostArray4D<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(IdefixHostArray4D<T>, _("IdefixHostArray4D<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if ( !convert and !py::array_t<T>::check_(src) )
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

        return true;
      }

      //Conversion part 2 (C++ -> Python)
      static py::handle cast(const IdefixHostArray4D<T>& src, py::return_value_policy policy, py::handle parent)
      {
        std::cout << "coucou2" << std::endl;
        std::vector<size_t> shape  (4);
        std::vector<size_t> strides(4);

        for ( int i = 0 ; i < 4 ; ++i ) {
          shape  [i] = src.extent  [i];
          strides[i] = src.strides[i]*sizeof(T);
        }

        //py::array a(std::move(shape), std::move(strides), src.data() );
        py::array_t<real, py::array::c_style> a({src.extent(0),src.extent(1),src.extent(2),src.extent(3)},src.data());
        return a.release();
      }
  };
}} // namespace pybind11::detail



/*
PYBIND11_PLUGIN(example) {
    py::module m("example", "Module description");
    m.def("func", &func, "Function description" );
    return m.ptr();
}*/
#endif // PYDEFIX_HPP_
