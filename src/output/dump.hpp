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
#include <vector>
#include <memory>
#include <type_traits>
#include <complex>
#include <cstdint>
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

// Underlying type is fixed explicitly: DataType is written as a raw int in
// dump files, so the on-disk format must not depend on the compiler default.
// The numeric values are part of the file format (pytools/dump_io.py relies on
// them): they MUST stay Double=0, Float=1, Int=2, Bool=3. New types are
// appended after to keep existing files readable.
enum class DataType : int {
    Double = 0,        /// <double>
    Float  = 1,        /// <float>
    Int    = 2,        /// <int>
    Bool   = 3,        /// <bool>
    ComplexFloat  = 4, /// <Kokkos::complex<float>>  (interleaved re,im floats)
    ComplexDouble = 5, /// <Kokkos::complex<double>> (interleaved re,im doubles)
    Int64         = 6, /// <int64_t>
    Unknown = -1       /// <Unknown type>
};
static_assert(static_cast<int>(DataType::Double) == 0, "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::Float) == 1, "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::Int) == 2, "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::Bool) == 3, "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::ComplexFloat) == 4,
              "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::ComplexDouble) == 5,
              "DataType value is on-disk format");
static_assert(static_cast<int>(DataType::Int64) == 6,
              "DataType value is on-disk format");

// Map C++ types to DataType
template<typename T> struct type_to_data_type {
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
template<> struct type_to_data_type<int64_t> {
    static constexpr DataType value = DataType::Int64;
};
// Kokkos::complex is layout-compatible with std::complex (two contiguous reals)
template<> struct type_to_data_type<std::complex<float>> {
    static constexpr DataType value = DataType::ComplexFloat;
};
template<> struct type_to_data_type<std::complex<double>> {
    static constexpr DataType value = DataType::ComplexDouble;
};
template<> struct type_to_data_type<Kokkos::complex<float>> {
    static constexpr DataType value = DataType::ComplexFloat;
};
template<> struct type_to_data_type<Kokkos::complex<double>> {
    static constexpr DataType value = DataType::ComplexDouble;
};

template<typename T>
constexpr DataType get_data_type() {
    return type_to_data_type<std::decay_t<T>>::value;
}

template<typename ArrayType>
constexpr DataType get_array_data_type() {
    // Kokkos::View<>::non_const_value_type
    using T = typename ArrayType::non_const_value_type;
    return get_data_type<T>();
}

// True when the array element type has a known on-disk representation
template<typename ArrayType>
inline constexpr bool scalar_or_known_type_v =
    get_array_data_type<ArrayType>() != DataType::Unknown;

// Size in bytes of a single element of the given type
constexpr int size_of(DataType type) {
    switch (type) {
        case DataType::Int:           return sizeof(int);
        case DataType::Float:         return sizeof(float);
        case DataType::Double:        return sizeof(double);
        case DataType::Bool:          return sizeof(bool);
        case DataType::ComplexFloat:  return 2 * sizeof(float);
        case DataType::ComplexDouble: return 2 * sizeof(double);
        case DataType::Int64:         return sizeof(int64_t);
        default:                      return 0;
    }
}

#ifdef WITH_MPI
// Matching MPI datatype (MPI_DATATYPE_NULL if unsupported)
inline MPI_Datatype to_mpi_type(DataType type) {
    switch (type) {
        case DataType::Int:           return MPI_INT;
        case DataType::Float:         return MPI_FLOAT;
        case DataType::Double:        return MPI_DOUBLE;
        case DataType::Bool:          return MPI_CXX_BOOL;
        case DataType::ComplexFloat:  return MPI_CXX_FLOAT_COMPLEX;
        case DataType::ComplexDouble: return MPI_CXX_DOUBLE_COMPLEX;
        case DataType::Int64:         return MPI_LONG_LONG;
        default:                      return MPI_DATATYPE_NULL;
    }
}
#endif

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

        // The storage semantics of a field follow from its rank: arrays with
        // rank <= 2 are replicated (identical on every rank), arrays of rank 3
        // or 4 are grid arrays and hence distributed. Raw pointers are scalar
        // and replicated. See agent/adr/0001.
        bool IsDistributed() const {
            return IsArray() && GetArrayDimension() >= 3;
        }

        // Shape of a replicated array (extent along each axis). Empty for
        // distributed fields and scalars.
        virtual std::vector<int> GetShape() const { return {}; }

        virtual void* GetRawPointer() const {
            IDEFIX_ERROR("GetRawPointer not implemented for this storage type");
            return nullptr;
        }

        // Generic host accessors used by the serial (replicated) path: a raw
        // pointer to host memory holding `Size()` elements of `GetDataType()`.
        virtual size_t Size() const { return 0; }
        // Copy this field's data into a host buffer of Size() elements
        virtual void CopyToHost(void* dest) const {
            IDEFIX_ERROR("CopyToHost not implemented for this storage type");
        }
        virtual void SyncFromHost(const void* source) const {
            IDEFIX_ERROR("SyncFromHost not implemented for this storage type");
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

        // Kokkos::View members used below: rank, non_const_value_type,
        // host_mirror_type, memory_space
        static constexpr int dimension = ArrayType::rank;
        using value_type = typename ArrayType::non_const_value_type;

        // Host mirror type with the same rank as the stored view
        using host_array_type = typename ArrayType::host_mirror_type;

        // True when the view's memory space is directly accessible from the
        // host (Kokkos::SpaceAccessibility, Kokkos_Concepts.hpp)
        static constexpr bool host_accessible =
            Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                       typename ArrayType::memory_space>::accessible;

        // 3D view of the stored data (directly for 3D arrays, else the var-th slice)
        auto GetSlice() const {
            if constexpr (dimension == 3) {
                return array;
            } else if constexpr (dimension == 4) {
                if (var < 0) IDEFIX_ERROR("var index required for 4D arrays");
                return Kokkos::subview(array, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
            } else {
                IDEFIX_ERROR("GetSlice is only defined for 3D/4D arrays");
            }
        }

     public:
        ArrayStorage(ArrayType& arr, int v, ArrayLocation loc, int dir)
            : array{arr}, var{v}, location{loc}, direction{dir} {}

        std::unique_ptr<StorageBase> clone() const override {
            return std::make_unique<ArrayStorage<ArrayType>>(array, var, location,
                                                             direction);
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
        int GetArrayDimension() const override { return dimension; }

        std::vector<int> GetShape() const override {
            std::vector<int> shape(dimension);
            for (int i = 0 ; i < dimension ; i++) shape[i] = array.extent(i);
            return shape;
        }

        // Total number of elements held by this rank (Kokkos::View::size())
        size_t Size() const override { return array.size(); }

        // Host copy of the whole array (rank-preserving), used by the serial path
        host_array_type GetHostArray() const {
            if constexpr (host_accessible) {
                return array;
            } else {
                auto hostArr = Kokkos::create_mirror(Kokkos::HostSpace(), array);
                Kokkos::deep_copy(hostArr, array);
                return hostArr;
            }
        }

        // Copy the array into an existing host buffer of Size() elements.
        // IdefixArrays are LayoutRight, hence contiguous in memory: the whole
        // array can be copied through a flat rank-1 host view (Kokkos::deep_copy)
        void CopyToHost(void* dest) const override {
            if constexpr (dimension == 1 || dimension == 2) {
                auto hostArr = GetHostArray();
                Kokkos::View<value_type*, Kokkos::HostSpace> dst(
                    static_cast<value_type*>(dest), Size());
                Kokkos::View<const value_type*, Kokkos::HostSpace> src(
                    hostArr.data(), Size());
                Kokkos::deep_copy(dst, src);
            } else {
                IDEFIX_ERROR("CopyToHost is only used for replicated 1D/2D arrays");
            }
        }

        void SyncFromHost(const void* source) const override {
            if constexpr (dimension == 1 || dimension == 2) {
                Kokkos::View<const value_type*, Kokkos::HostSpace> src(
                    static_cast<const value_type*>(source), Size());
                if constexpr (host_accessible) {
                    Kokkos::View<value_type*, Kokkos::HostSpace> dst(array.data(), Size());
                    Kokkos::deep_copy(dst, src);
                } else {
                    auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), array);
                    Kokkos::View<value_type*, Kokkos::HostSpace> dst(mirror.data(), Size());
                    Kokkos::deep_copy(dst, src);
                    Kokkos::deep_copy(array, mirror);
                }
            } else {
                IDEFIX_ERROR("SyncFromHost is only used for replicated 1D/2D arrays");
            }
        }

        void* GetRawPointer() const override {
            IDEFIX_ERROR("GetRawPointer not supported for array storage");
            return nullptr;
        }

        IdefixHostArray3D<real> GetHostField() const override {
            if constexpr (dimension == 3 || dimension == 4) {
                if constexpr (std::is_same_v<value_type, real>) {
                    if constexpr (host_accessible) {
                        return GetSlice();
                    } else {
                        auto hostArr = Kokkos::create_mirror(Kokkos::HostSpace(), GetSlice());
                        Kokkos::deep_copy(hostArr, GetSlice());
                        return hostArr;
                    }
                } else {
                    IDEFIX_ERROR("GetHostField requires real-valued arrays");
                    return IdefixHostArray3D<real>();
                }
            } else {
                IDEFIX_ERROR("GetHostField is only defined for 3D/4D arrays");
                return IdefixHostArray3D<real>();
            }
        }

        void SyncFrom(const IdefixHostArray3D<real>& source) const override {
            if constexpr (dimension != 3 && dimension != 4) {
                IDEFIX_ERROR("Only 3D and 4D arrays supported");
            } else if constexpr (host_accessible) {
                Kokkos::deep_copy(GetSlice(), source);
            } else {
                auto mirror = Kokkos::create_mirror(Kokkos::HostSpace(), GetSlice());
                Kokkos::deep_copy(mirror, source);
                Kokkos::deep_copy(GetSlice(), mirror);
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

        size_t Size() const override { return static_cast<size_t>(size); }

        // Raw storage already lives in host memory
        void CopyToHost(void* dest) const override {
            T* dst = static_cast<T*>(dest);
            for (int n = 0 ; n < size ; n++) dst[n] = data[n];
        }

        void SyncFromHost(const void* source) const override {
            const T* src = static_cast<const T*>(source);
            for (int n = 0 ; n < size ; n++) data[n] = src[n];
        }

        void* GetRawPointer() const override {
            return static_cast<void*>(data);
        }

        IdefixHostArray3D<real> GetHostField() const override {
            IDEFIX_ERROR("GetHostField not supported for raw storage - use GetRawPointer instead");
            return IdefixHostArray3D<real>();
        }

        void SyncFrom(const IdefixHostArray3D<real>& source) const override {
            // Raw storage already aliases the caller's data: nothing to copy back
        }
    };

 public:
    // =============================================================================
    // Constructors
    // =============================================================================

    // IdefixArray storage: the storage semantics follow from the rank, so no
    // option is offered to the caller (see agent/adr/0001).
    template<typename ArrayType>
    DumpField(ArrayType& arr, int varnum = -1,
              ArrayLocation loc = ArrayLocation::Center, int dir = -1)
        : impl{std::make_unique<ArrayStorage<ArrayType>>(arr, varnum, loc, dir)} {}

    // Scalar (raw pointer) storage: any type with a known DataType, including
    // std::complex<float/double> and Kokkos::complex<float/double>.
    template<typename T>
    explicit DumpField(T* data, int size = 1)
        : impl{std::make_unique<RawStorage<T>>(data, size)} {
        static_assert(dump_io::get_data_type<T>() != dump_io::DataType::Unknown,
                      "Unsupported scalar type for dump I/O");
    }

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

    template<typename T>
    T* GetRawPointer() const {
        return static_cast<T*>(impl->GetRawPointer());
    }

    void SyncFrom(const IdefixHostArray3D<real>& source) const {
        impl->SyncFrom(source);
    }

    // Serial (replicated) path: copy this field to/from a host buffer of
    // Size() elements, laid out contiguously in the field's DataType.
    void CopyToHost(void* dest) const { impl->CopyToHost(dest); }
    void SyncFromHost(const void* source) const { impl->SyncFromHost(source); }

    // =============================================================================
    // Accessors
    // =============================================================================

    ArrayLocation GetLocation() const { return impl->GetLocation(); }
    int GetDirection() const { return impl->GetDirection(); }
    int GetSize() const { return impl->GetSize(); }
    dump_io::DataType GetDataType() const { return impl->GetDataType(); }
    bool IsArray() const { return impl->IsArray(); }
    int GetVarNum() const { return impl->GetVarNum(); }
    // 3D/4D grid arrays are distributed across ranks; 1D/2D arrays and scalars
    // are replicated (identical on every rank).
    bool IsDistributed() const { return impl->IsDistributed(); }
    std::vector<int> GetShape() const { return impl->GetShape(); }
    size_t Size() const { return impl->Size(); }
    int GetArrayDimension() const { return impl->GetArrayDimension(); }

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
    // Storage is decided by the array rank, with no override: 3D/4D grid arrays
    // are distributed across ranks, 1D/2D arrays are replicated (identical on
    // every rank). See agent/adr/0001.
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
    std::vector<char> serialBuffer;         // Host buffer for replicated fields

    std::map<std::string, DumpField> dumpFieldMap;

    // Copy a replicated field to a host buffer of bytes, resizing it as needed
    void* GetSerialBuffer(DumpField const&);

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
    static_assert(Kokkos::is_view_v<ArrayType>,
                 "RegisterVariable requires an IdefixArray type");
    static_assert(ArrayType::rank >= 1 && ArrayType::rank <= 4,
                 "Only 1D, 2D, 3D and 4D arrays can be dumped");
    // Distributed (3D/4D) grid fields are stored as real-precision sub-arrays.
    // Replicated (1D/2D) fields may hold any type supported by dump_io::DataType.
    static_assert(ArrayType::rank <= 2 ||
                  std::is_same_v<typename ArrayType::non_const_value_type, real>,
                  "Only real-valued arrays can be dumped as distributed fields");
    static_assert(dump_io::scalar_or_known_type_v<ArrayType>,
                  "Unsupported array element type for dump I/O");
    // A 4D array stores several variables: varnum designates which one to dump.
    // Other ranks have a single variable and ignore varnum.
    if (ArrayType::rank != 4 && varnum >= 0) {
        IDEFIX_ERROR("varnum is only meaningful for 4D arrays");
    }

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
