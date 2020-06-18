# TODO: One could speed up compilation by compiling only the module compatible with HD/MHD
ifdef MHD
INCLUDE_DIR += -I$(SRC)/hydro/MHDsolvers
VPATH += $(SRC)/hydro/MHDsolvers
HEADERS += solversMHD.hpp
OBJ +=  hlldMHD.o hllMHD.o roeMHD.o tvdlfMHD.o
else
INCLUDE_DIR += -I$(SRC)/hydro/HDsolvers
VPATH += $(SRC)/hydro/HDsolvers
HEADERS += solversHD.hpp
OBJ += hllcHD.o hllHD.o roeHD.o tvdlfHD.o
endif
