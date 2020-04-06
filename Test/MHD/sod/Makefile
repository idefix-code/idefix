
SRC 		 = $(IDEFIX_DIR)/src
INCLUDE_DIRS = -I, -I. -I$(SRC) -I$(SRC)/MHD -I$(SRC)/Test
VPATH 		 = ./:$(SRC):$(SRC)/MHD:$(SRC)/Test

KOKKOS_PATH = ${IDEFIX_DIR}/src/kokkos
KOKKOS_DEVICES = "Serial"
KOKKOS_ARCH = "Pascal61"

# include local rules which are not hosted on git
# Can Makefile.local can redefine the make rules
-include Makefile.local

EXE_NAME = "idefix"




HEADERS = definitions.hpp error.hpp arrays.hpp grid.hpp setup.hpp gridHost.hpp idefix.hpp dataBlock.hpp dataBlockHost.hpp mod_defs.hpp input.hpp kokkos_types.h loop.hpp real_types.h timeIntegrator.hpp test.hpp outputVtk.hpp physics.hpp
OBJ = grid.o gridHost.o input.o main.o dataBlock.o dataBlockHost.o setup.o test.o timeIntegrator.o outputVtk.o physics.o error.o

#VPATH="src/"

default: build


ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
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

# Create gitversion if Idefix is in a git repo
ifneq ("$(wildcard $(IDEFIX_DIR)/.git/HEAD))","")
GIT_DEPENDENCE = $(IDEFIX_DIR)/.git/HEAD $(IDEFIX_DIR)/.git/index
else
GIT_DEPENDENCE =
endif

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
$(SRC)/gitversion.h: $(GIT_DEPENDENCE)
	echo "#define GITVERSION \"$(shell cd $(IDEFIX_DIR); git describe --tags --always)\"" > $@

#specific dependence of input.o
input.o: $(SRC)/gitversion.h $(HEADERS) input.cpp

test: $(EXE)
	./$(EXE)
