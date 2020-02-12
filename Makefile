

KOKKOS_PATH = ${HOME}/Kokkos/kokkos
KOKKOS_DEVICES = "OpenMP"

# include local rules which are not hosted on git
# Can Makefile.local can redefine the make rules
-include Makefile.local

EXE_NAME = "idefix"


IDEFIX_DIR 	 = ./
SRC 		 = $(IDEFIX_DIR)/src
INCLUDE_DIRS = -I, -I. -I$(SRC) -I$(SRC)/HD -I$(SRC)/Test
VPATH 		 = ./:$(SRC):$(SRC)/HD:$(SRC)/Test

HEADERS = arrays.hpp globals.hpp grid.hpp gridHost.hpp idefix.hpp dataBlock.hpp dataBlockHost.hpp mod_defs.hpp input.hpp kokkos_types.h loop.hpp real_types.h timeIntegrator.hpp test.hpp
OBJ = globals.o grid.o gridHost.o input.o main.o dataBlock.o dataBlockHost.o test.o timeIntegrator.o

#VPATH="src/"

default: build


ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "Pascal61"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "BDW"
endif

CXXFLAGS = -O3 -g
LINK = ${CXX}
LINKFLAGS = 

DEPFLAGS = -M

LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) $(INCLUDE_DIRS) -c $<

#dependance on headers
$(OBJ): $(HEADERS)

# in order to have access to current git version
$(SRC)/gitversion.h: .git/HEAD .git/index
	echo "#define GITVERSION \"$(shell git describe --tags --always)\"" > $@

#specific dependence of input.o
input.o: $(SRC)/gitversion.h

test: $(EXE)
	./$(EXE)
