// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "idefix.hpp"
#include "readCSV.hpp"

ReadCSV::ReadCSV(std::string filename, char delimiter) {
  idfx::pushRegion("ReadCSV::ReadCSV");
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
        std::vector<real> dataLine;
        dataLine.clear();
        // make the line a string stream, and get all of the values separated by a delimiter
        std::stringstream str(line);
        std::string valueString;
        while(std::getline(str, valueString, delimiter)) {
          if(firstLine && firstColumn) {
            // We skip this
            firstColumn = false;
            continue;
          }
          real value;
          try {
            value = std::stod(valueString);
          } catch(const std::exception& e) {
            std::stringstream errmsg;
            errmsg << e.what() << std::endl
                   << "ReadCSV: Error while parsing " << filename  << ", \"" << valueString
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
            IDEFIX_ERROR("ReadCSV: The number of columns in the input CSV file should be constant");
          }
          dataVector.push_back(dataLine);
          firstLine = false;
        }
      }
      file.close();
      // End of file reached
    } else {
      std::stringstream errmsg;
      errmsg << "ReadCSV: Unable to open file " << filename << std::endl;
      IDEFIX_ERROR(errmsg);
    }

    size[0] = xVector.size();
    size[1] = yVector.size();
  }

  #ifdef WITH_MPI
    // Share the size of the arrays
    MPI_Bcast(size, 2, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  //Allocate arrays so that the data fits in it
  this->x = IdefixArray1D<real> ("CSV_x", size[0]);
  this->y = IdefixArray1D<real> ("CSV_y", size[1]);
  this->data =  IdefixArray2D<real> ("CSV_data", size[0], size[1]);

  IdefixArray1D<real>::HostMirror xHost = Kokkos::create_mirror_view(this->x);
  IdefixArray1D<real>::HostMirror yHost = Kokkos::create_mirror_view(this->y);
  IdefixArray2D<real>::HostMirror dataHost = Kokkos::create_mirror_view(this->data);

  // Fill the arrays with the std::vector content
  if(idfx::prank == 0) {
    for(int i = 0 ; i < xVector.size(); i++) {
      xHost(i) = xVector[i];
    }
    for(int i = 0 ; i < yVector.size(); i++) {
      yHost(i) = yVector[i];
    }
    for(int i = 0 ; i < dataVector.size(); i++) {
      auto line = dataVector[i];
      for(int j = 0 ; j < line.size(); j++) {
        dataHost(i,j) = line[j];
      }
    }
  }

  #ifdef WITH_MPI
    // Share with the others
    MPI_Bcast(xHost.data(), size[0], realMPI, 0, MPI_COMM_WORLD);
    MPI_Bcast(yHost.data(), size[1], realMPI, 0, MPI_COMM_WORLD);
    MPI_Bcast(dataHost.data(), size[0]*size[1], realMPI, 0, MPI_COMM_WORLD);
  #endif

  // Copy to target
  Kokkos::deep_copy(this->x ,xHost);
  Kokkos::deep_copy(this->y, yHost);
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

KOKKOS_FUNCTION real K_ReadCSVFetch(const real x, const real y,
                                    const IdefixArray1D<real> &xin,
                                    const IdefixArray1D<real> &yin,
                                    const IdefixArray2D<real> &datain) {
  int nx = xin.extent(0);
  int ny = yin.extent(0);

  int endx = nx-1;
  int endy = ny-1;
  if(x < xin(0) || x >xin(endx)) {
    return(NAN);
  }
  if(y < yin(0) || y >yin(endy)) {
    return(NAN);
  }

  int i = static_cast<int> ( (x - xin(0)) / (xin(endx) - xin(0)) * nx );
  int j = static_cast<int> ( (y - yin(0)) / (yin(endx) - yin(0)) * ny );

  if(xin(i) > x || xin(i+1) < x) { // x-points are not evenly distributed
    i = 0;
    while(xin(i) > x) {
      i++;
    }
  }
  if(yin(j) > y || yin(j+1) < y) { // y-points are not evenly distributed
    j = 0;
    while(yin(j) > y) {
      j++;
    }
  }

  real delta_x = (x - xin(i)) / (xin(i+1) - xin(i));
  real delta_y = (y - yin(j)) / (yin(j+1) - yin(j));

  // Rectilinear interpolation, in principle only valid for regularly-spaced datapoints.
  real value =   (1.0-delta_x) * (1.0-delta_y) * datain(i  ,j  )
                       + (1.0-delta_x) * (delta_y)     * datain(i  ,j+1)
                       + (delta_x)     * (1.0-delta_y) * datain(i+1,j  )
                       + (delta_x)     * (delta_y)     * datain(i+1,j+1);



  return(value);
}
