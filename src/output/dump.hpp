// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/output/dump.hpp
//
// Last modified : 11/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2026)
// - Clément Robert <clement.robert@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2020 - 2024)
// - Gaylor Wafflard <gaylor.wafflard@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2024)
// - Alankar Dutta <dutta.alankar@gmail.com> (MPA Garching - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef OUTPUT_DUMP_HPP_
#define OUTPUT_DUMP_HPP_
#include <string>
#include <map>
#include <array>
#include <memory>
#include <type_traits>
#if __has_include(<filesystem>)
  #include <filesystem> // NOLINT [build/c++17]
  namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#else
  #error "Missing the <filesystem> header."
#endif
#include "idefix.hpp"
#include "input.hpp"
#include "dataBlock.hpp"
#include "arrays.hpp"

// Define data descriptor used for distributed I/O when MPI is enabled
#ifdef WITH_MPI
  using IdfxDataDescriptor = MPI_Datatype;
#else
  using IdfxDataDescriptor = int;   // This is actually not used
#endif

// Forward class declaration
//class Vtk;
class Output;
class DataBlock;

// =====================================================================================
// I/O-specific DataType enumeration
// =====================================================================================

namespace dump_io {

enum class DataType {
    Int,           /// <int>
    Float,         /// <float>
    Double,        /// <double>
    Bool,          /// <bool>
    ComplexFloat,  /// <std::complex<float>>
    ComplexDouble, /// <std::complex<double>>
    Unknown        /// <Unknown type>
};

// Map C++ types to DataType
template<typename T, typename = void> struct type_to_data_type {
    static constexpr DataType value = DataType::Unknown;
};
template<> struct type_to_data_type<int> {
    static constexpr DataType value = DataType::Int;
};
template<> struct type_to_data_type<float> {
    static constexpr DataType value = DataType::Float;
};
template<> struct type_to_data_type<double> {
    static constexpr DataType value = DataType::Double;
};
template<> struct type_to_data_type<bool> {
    static constexpr DataType value = DataType::Bool;
};
template<> struct type_to_data_type<std::complex<float>> {
    static constexpr DataType value = DataType::ComplexFloat;
};
template<> struct type_to_data_type<std::complex<double>> {
    static constexpr DataType value = DataType::ComplexDouble;
};

template<typename T>
constexpr DataType get_data_type() {
    return type_to_data_type<std::decay_t<T>>::value;
}

template<typename ArrayType>
constexpr DataType get_array_data_type() {
    using T = idefix::value_type_t<ArrayType>;
    return get_data_type<T>();
}

} // namespace dump_io

// =====================================================================================
// GridBox structure
// =====================================================================================

struct GridBox {
    std::array<int,3> start;
    std::array<int,3> size;
    std::array<int,3> sizeGlob;
};

// =====================================================================================
// Storage classes for DumpField (defined first so they can be used in DumpField)
// =====================================================================================

class DumpField {
 public:
    enum class ArrayLocation { Center, Face, Edge };

 private:
    // =============================================================================
    // Base class for type-erased storage
    // =============================================================================

    class StorageBase {
     public:
        virtual ~StorageBase() = default;
        virtual std::unique_ptr<StorageBase> clone() const = 0;

        virtual ArrayLocation GetLocation() const = 0;
        virtual int GetDirection() const = 0;
        virtual int GetSize() const = 0;
        virtual dump_io::DataType GetDataType() const = 0;
        virtual bool IsArray() const = 0;
        virtual int GetVarNum() const = 0;
        virtual int GetArrayDimension() const = 0;
        virtual bool IsHost() const = 0;
        virtual void SyncToHost() const = 0;

        virtual void* GetRawPointer() const {
            IDEFIX_ERROR("GetRawPointer not implemented for this storage type");
            return nullptr;
        }

        virtual IdefixHostArray3D<real> GetHostField() const {
            IDEFIX_ERROR("GetHostField not implemented for this storage type");
            return IdefixHostArray3D<real>();
        }

        virtual void SyncFrom(const IdefixHostArray3D<real>& source) const {
            IDEFIX_ERROR("SyncFrom not implemented for this storage type");
        }
    };

    // =============================================================================
    // Storage implementation for IdefixArrays
    // =============================================================================

    template<typename ArrayType>
    class ArrayStorage : public StorageBase {
        ArrayType& array;
        int var;
        ArrayLocation location;
        int direction;

     public:
        ArrayStorage(ArrayType& arr, int v, ArrayLocation loc, int dir)
            : array{arr}, var{v}, location{loc}, direction{dir} {}

        std::unique_ptr<StorageBase> clone() const override {
            return std::make_unique<ArrayStorage<ArrayType>>(array, var, location, direction);
        }

        // Accessors
        ArrayLocation GetLocation() const override { return location; }
        int GetDirection() const override { return direction; }
        int GetSize() const override { return -1; }
        dump_io::DataType GetDataType() const override {
            return dump_io::get_array_data_type<ArrayType>();
        }
        bool IsArray() const override { return true; }
        int GetVarNum() const override { return var; }
        int GetArrayDimension() const override {
            return idefix::view_dimension_v<ArrayType>;
        }
        bool IsHost() const override {
            return idefix::is_host_space_v<ArrayType>;
        }

        void* GetRawPointer() const override {
            IDEFIX_ERROR("GetRawPointer not supported for array storage");
            return nullptr;
        }

        IdefixHostArray3D<real> GetHostField() const override {
            using T = idefix::value_type_t<ArrayType>;
            static_assert(std::is_same_v<T, real>, "GetHostField requires real type for now");

            if constexpr (idefix::view_dimension_v<ArrayType> == 3) {
                if constexpr (idefix::is_host_space_v<ArrayType>) {
                    return array;
                } else {
                    IdefixHostArray3D<real> hostArr =
                        Kokkos::create_mirror(Kokkos::HostSpace(), array);
                    Kokkos::deep_copy(hostArr, array);
                    return hostArr;
                }
            }
            else if constexpr (idefix::view_dimension_v<ArrayType> == 4) {
                if (var < 0) {
                    IDEFIX_ERROR("var index required for 4D arrays");
                }
                auto slice = Kokkos::subview(array, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
                if constexpr (idefix::is_host_space_v<ArrayType>) {
                    return slice;
                } else {
                    IdefixHostArray3D<real> hostArr =
                        Kokkos::create_mirror(Kokkos::HostSpace(), slice);
                    Kokkos::deep_copy(hostArr, slice);
                    return hostArr;
                }
            }
            else {
                static_assert(idefix::view_dimension_v<ArrayType> == 3 ||
                             idefix::view_dimension_v<ArrayType> == 4,
                             "Only 3D and 4D arrays supported");
            }
        }

        template<typename SourceArrayType>
        void SyncFrom(const SourceArrayType& source) const {
            // Copy from source array to our stored array
            if constexpr (idefix::view_dimension_v<ArrayType> == 3) {
                if constexpr (idefix::view_dimension_v<SourceArrayType> == 3) {
                    Kokkos::deep_copy(array, source);
                }
            }
            else if constexpr (idefix::view_dimension_v<ArrayType> == 4) {
                if constexpr (idefix::view_dimension_v<SourceArrayType> == 3) {
                    // Source is 3D, we need to copy to the var-th slice of our 4D array
                    auto destSlice =
                        Kokkos::subview(array, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
                    Kokkos::deep_copy(destSlice, source);
                }
            }
        }

        void SyncToHost() const override {
            if constexpr (!idefix::is_host_space_v<ArrayType>) {
                if constexpr (idefix::view_dimension_v<ArrayType> == 4) {
                    if (var >= 0) {
                        auto slice =
                            Kokkos::subview(array, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
                        auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), slice);
                        Kokkos::deep_copy(mirror, slice);
                    }
                } else {
                    auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), array);
                    Kokkos::deep_copy(mirror, array);
                }
            }
        }

        void SyncFrom(const IdefixHostArray3D<real>& source) const override {
            // Copy from source array to our stored array
            if constexpr (idefix::view_dimension_v<ArrayType> == 3) {
                if constexpr (idefix::is_host_space_v<ArrayType>) {
                    Kokkos::deep_copy(array, source);
                } else {
                    // For device arrays, we need to copy from host to device
                    auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), array);
                    Kokkos::deep_copy(mirror, source);
                    Kokkos::deep_copy(array, mirror);
                }
            }
            else if constexpr (idefix::view_dimension_v<ArrayType> == 4) {
                if (var >= 0) {
                    auto destSlice =
                        Kokkos::subview(array, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
                    if constexpr (idefix::is_host_space_v<ArrayType>) {
                        Kokkos::deep_copy(destSlice, source);
                    } else {
                        // For device 4D arrays, copy to the var-th slice
                        auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), destSlice);
                        Kokkos::deep_copy(mirror, source);
                        Kokkos::deep_copy(destSlice, mirror);
                    }
                }
            }
        }
    };

    // =============================================================================
    // Storage implementation for raw pointers
    // =============================================================================

    template<typename T>
    class RawStorage : public StorageBase {
        T* data;
        int size;

     public:
        RawStorage(T* d, int s) : data{d}, size{s} {}

        std::unique_ptr<StorageBase> clone() const override {
            return std::make_unique<RawStorage<T>>(data, size);
        }

        ArrayLocation GetLocation() const override { return ArrayLocation::Center; }
        int GetDirection() const override { return -1; }
        int GetSize() const override { return size; }
        dump_io::DataType GetDataType() const override {
            return dump_io::get_data_type<T>();
        }
        bool IsArray() const override { return false; }
        int GetVarNum() const override { return -1; }
        int GetArrayDimension() const override { return 0; }
        bool IsHost() const override { return true; }

        void* GetRawPointer() const override {
            return static_cast<void*>(data);
        }

        IdefixHostArray3D<real> GetHostField() const override {
            IDEFIX_ERROR("GetHostField not supported for raw storage - use GetRawPointer instead");
            return IdefixHostArray3D<real>();
        }

        void SyncToHost() const override {
            // Nothing to do - already on host
        }

        void SyncFrom(const IdefixHostArray3D<real>& source) const override {
            // For raw pointers, we copy directly from the source
            // This is a no-op for raw storage as the data is already in the right place
        }
    };

 public:
    // =============================================================================
    // Constructors
    // =============================================================================

    template<typename ArrayType>
    DumpField(ArrayType& arr, int varnum = -1,
              ArrayLocation loc = ArrayLocation::Center, int dir = -1)
        : impl{std::make_unique<ArrayStorage<ArrayType>>(arr, varnum, loc, dir)} {}

    template<typename T>
    explicit DumpField(T* data, int size = 1)
        : impl{std::make_unique<RawStorage<T>>(data, size)} {}

    // =============================================================================
    // Copy and move semantics
    // =============================================================================

    DumpField(const DumpField& other) : impl{other.impl->clone()} {}
    DumpField(DumpField&& other) noexcept = default;
    DumpField& operator=(const DumpField& other) {
        if (this != &other) impl = other.impl->clone();
        return *this;
    }
    DumpField& operator=(DumpField&& other) noexcept = default;
    ~DumpField() = default;

    // =============================================================================
    // Access methods
    // =============================================================================

    // Default version returns real
    IdefixHostArray3D<real> GetHostField() const {
        return impl->GetHostField();
    }

    // Template version for backward compatibility
    template<typename T>
    T GetHostField() const {
        if constexpr (std::is_same_v<T, IdefixHostArray3D<real>>) {
            return impl->GetHostField();
        } else if constexpr (std::is_pointer_v<T>) {
            // For pointer types, return the raw pointer cast to the requested type
            return static_cast<T>(impl->GetRawPointer());
        } else {
            // For non-pointer types, dereference and return
            return *static_cast<T*>(impl->GetRawPointer());
        }
    }

    void* GetRawPointer() const {
        return impl->GetRawPointer();
    }

    template<typename T>
    T* GetRawPointer() const {
        return static_cast<T*>(impl->GetRawPointer());
    }

    void SyncToHost() const {
        impl->SyncToHost();
    }

    void SyncFrom(const IdefixHostArray3D<real>& source) const {
        impl->SyncFrom(source);
    }

    // =============================================================================
    // Accessors
    // =============================================================================

    ArrayLocation GetLocation() const { return impl->GetLocation(); }
    int GetDirection() const { return impl->GetDirection(); }
    int GetSize() const { return impl->GetSize(); }
    dump_io::DataType GetDataType() const { return impl->GetDataType(); }
    bool IsArray() const { return impl->IsArray(); }
    int GetVarNum() const { return impl->GetVarNum(); }
    int GetArrayDimension() const { return impl->GetArrayDimension(); }
    bool IsHost() const { return impl->IsHost(); }

 private:
    // =============================================================================
    // Member variable - type-erased storage
    // =============================================================================
    std::unique_ptr<StorageBase> impl;
};

// =====================================================================================
// Dump class declaration
// =====================================================================================

class Dump {
    friend class DumpImage; // Allow dumpimag to have access to dump API
    friend class Output;

 public:
    explicit Dump(Input &, DataBlock *);               // Create Dump Object
    explicit Dump(DataBlock *);               // Create a dump object independent of input
    ~Dump();

    // Create a Dump file from the current state of the code
    int Write(Output&);
    // Read and load a dump file as current state of the code
    bool Read(Output&, int);

    // Register IdefixArrays - generic versions
    template<typename ArrayType>
    void RegisterVariable(ArrayType& in, std::string name,
                         int varnum = -1, int dir = -1,
                         DumpField::ArrayLocation loc = DumpField::ArrayLocation::Center);

    // Register any other fundamental type
    template<typename T>
    void RegisterVariable(T* in, std::string name, int size = 1);

 private:
    void Init(DataBlock*);
    DataBlock *data;
    int dumpFileNumber;
    int geometry{GEOMETRY};
    int periodicity[3];

    real *scrch;                            // Scratch array in host space

    std::map<std::string, DumpField> dumpFieldMap;

    // Timer
    Kokkos::Timer timer;

    // File offset
#ifdef WITH_MPI
    MPI_Offset offset;
#endif
    // These descriptors are only useful with MPI
    IdfxDataDescriptor descCR;   // Descriptor for cell-centered fields (Read)
    IdfxDataDescriptor descCW;   // Descriptor for cell-centered fields (Write)
    IdfxDataDescriptor descSR[3]; // Descriptor for face-centered fields (Read)
    IdfxDataDescriptor descSW[3]; // Descriptor for face-centered fields (Write)
    IdfxDataDescriptor descER[3]; // Descriptor for edge-centered fields (Read)
    IdfxDataDescriptor descEW[3]; // Descriptor for edge-centered fields (Write)

    void WriteString(IdfxFileHandler, char *, int);
    void WriteSerial(IdfxFileHandler, int, int *, dump_io::DataType, char*, void*);
    void WriteDistributed(IdfxFileHandler, int, int*, int*, char*, IdfxDataDescriptor&, real*);
    void ReadNextFieldProperties(IdfxFileHandler, int&, int*, dump_io::DataType&, std::string&);
    void ReadSerial(IdfxFileHandler, int, int*, dump_io::DataType, void*);
    void ReadDistributed(IdfxFileHandler, int, int*, int*, IdfxDataDescriptor&, void*);
    void Skip(IdfxFileHandler, int, int *, dump_io::DataType);
    int GetLastDumpInDirectory(fs::path &);
    void CreateMPIDataType(GridBox, bool);

    fs::path outputDirectory;
};

// =====================================================================================
// RegisterVariable template implementations
// =====================================================================================

template<typename ArrayType>
void Dump::RegisterVariable(ArrayType& in, std::string name,
                           int varnum, int dir,
                           DumpField::ArrayLocation loc) {
    static_assert(idefix::is_idefix_array_v<ArrayType>,
                 "RegisterVariable requires an IdefixArray type");

    if (dumpFieldMap.count(name) > 0) {
        std::stringstream msg;
        msg << "Error while registering variable in Dump I/O: " << std::endl
            << name << " has already been registered." << std::endl;
        IDEFIX_ERROR(msg);
    }
    dumpFieldMap.emplace(name, DumpField(in, varnum, loc, dir));
}

template<typename T>
void Dump::RegisterVariable(T* in, std::string name, int size) {
    if (dumpFieldMap.count(name) > 0) {
        std::stringstream msg;
        msg << "Error while registering variable in Dump I/O: " << std::endl
            << name << " has already been registered." << std::endl;
        IDEFIX_ERROR(msg);
    }
    dumpFieldMap.emplace(name, DumpField(in, size));
}

#endif // OUTPUT_DUMP_HPP_
