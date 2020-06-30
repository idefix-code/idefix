
SRC 		 = $(IDEFIX_DIR)/src
INCLUDE_DIRS = -I, -I. -I$(SRC) -I$(SRC)/hydro -I$(SRC)/Test  -I$(SRC)/hydro/HDsolvers
VPATH 		 = ./:$(SRC):$(SRC)/hydro:$(SRC)/Test:$(SRC)/hydro/HDsolvers

KOKKOS_PATH = ${IDEFIX_DIR}/src/kokkos
KOKKOS_DEVICES = "Serial"
KOKKOS_ARCH = "BDW"

# include local rules which are not hosted on git
# Can Makefile.local can redefine the make rules
-include Makefile.local

EXE = "idefix"

HEADERS = definitions.hpp error.hpp arrays.hpp grid.hpp setup.hpp gridHost.hpp idefix.hpp dataBlock.hpp dataBlockHost.hpp global.hpp input.hpp kokkos_types.h loop.hpp real_types.h timeIntegrator.hpp test.hpp outputVtk.hpp hydro.hpp
HEADERS +=  solversHD.hpp
OBJ = grid.o gridHost.o input.o main.o dataBlock.o global.o dataBlockHost.o setup.o test.o timeIntegrator.o outputVtk.o hydro.o error.o
OBJ +=  hllcHD.o hllHD.o roeHD.o tvdlfHD.o

#VPATH="src/"

default: build

CXX = mpicxx
CXXFLAGS =  -O3 -D WITH_MPI -D MHD=NO



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
