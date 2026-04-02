/**
 * @file VTI.cpp
 * @brief DAMASK VTI reader with beast-based XML parser
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#ifdef BOOST

#include <cstdint>
#include <cstddef>
#include <format> // IWYU pragma: keep
#include <string>
#include <string_view>
#include <vector>
#include "ISO_Fortran_binding.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include "../IO.h"

// Guideline Support Library is used when pointers own memory and need to be manually freed.
// https://clang.llvm.org/extra/clang-tidy/checks/cppcoreguidelines/owning-memory.html
namespace gsl {
  template<typename T>
  using owner = T;
}

namespace pt = boost::property_tree;

struct DecodedBuffer {
    std::string vtk_type;
    std::vector<uint8_t> raw_bytes;
};

/**
 * @brief Holds the parsed VTI property tree.
 */
class VTI {
public:
  /**
   * @brief Construct by reading and parsing a VTI file from disk.
   * @param file_path Null-terminated path to a .vti file.
   * @throws std::runtime_error on I/O or parse/validation errors.
   */
  explicit VTI(const IO& io) : io(io) {}
  VTI(const char* file_path, const IO& io);

  /**
   * @brief Read a little-endian word of size @p n_bytes_per_word from @p p.
   *
   * @param p           Pointer to the first byte of the word.
   * @param n_bytes_per_word  Width of the word (4 or 8).
   * @return            Zero-extended 64-bit unsigned integer result.
   */
  static uint64_t read_word(const uint8_t* p, std::size_t n_bytes_per_word);

  /**
   * @brief Inflate and concatenate all compressed blocks in a VTI DataArray.
   *
   * @param b64_string  Raw Base-64 string
   * @param n_bytes_per_word Word size (4 or 8).
   * @return            Vector with decoded uncompressed bytes
   */
  std::vector<uint8_t> decode_compressed_VTI(const std::string& b64_string, std::size_t n_bytes_per_word);

  /**
   * @brief Decode an uncompressed VTI Dataarray.
   *
   * @param b64_string  Raw Base-64 string.
   * @param n_bytes_per_word Word size (4 or 8).
   * @return            Vector with decoded bytes.
   */
  std::vector<uint8_t> decode_uncompressed_VTI(const std::string& b64_string, std::size_t n_bytes_per_word);

  /**
   * @brief Read an integer DataArray into a Fortran pointer descriptor.
   *
   * @param name  Name of target attribute.
   * @param desc  Pre-allocated descriptor to be filled by \c CFI_allocate.
   */
  void read_dataset_int(const char* name, CFI_cdesc_t* desc);

  /**
   * @brief Read a floating-point DataArray into a Fortran pointer descriptor.
   *
   * @param name  Name of target attribute.
   * @param desc  Pre-allocated descriptor to be filled by \c CFI_allocate.
   */
  void read_dataset_real(const char* name, CFI_cdesc_t* desc);

  /**
   * @brief Extract grid size, physical extent and origin from a VTI file.
   *
   * @param cells_ptr     Number of cells along x/y/z.
   * @param geom_size_ptr  Physical side lengths.
   * @param origin_ptr    Origin coordinates.
   * @param labels_desc   Optional labels descriptor (may be nullptr).
   */
  void read_geometry(int* cells_ptr,
                     double* geom_size_ptr,
                     // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                     double* origin_ptr,
                     CFI_cdesc_t* labels_desc);

  /**
   * @brief Locate a VTK DataArray inside the class VTKFile buffer and return its bytes.
   *
   * @param array_name Name of the target DataArray.
   * @return           Struct with vtk datatype and the raw decoded bytes.
   */
  DecodedBuffer parse_cell_data_array(const char* array_name);
  pt::ptree vti_tree;
  std::string file_path;
  IO io;

private:
  /**
   * @brief Fetch an XML attribute, return an empty string if it doesn't exist.
   *
   * @param n    XML node.
   * @param key  Attribute name.
   * @return     A string with the attribute if it exists, otherwise an empty one.
   */
  std::string get_attr(const pt::ptree& n, const char* key) const;

  /**
   * @brief Decodes a Base-64 string using Boost.Beast.
   *
   * @param b64 ASCII string with valid Base-64 characters
   * @return Vector with decoded bytes
   */
  std::vector<std::uint8_t> decode_b64(std::string_view b64);

  /**
   * @brief Interpret a raw byte vector as a span of type T.
   *
   * @tparam T   Destination type.
   * @param raw  Vector of raw bytes
   * @return     Span of type T pointing into @p raw.
   */
  template<class T>
  std::vector<T> view(const std::vector<uint8_t>& raw) const;

  /**
   * @brief Allocate a Fortran array and convert the VTK byte stream into it.
   *
   * @tparam T   Target element type (int or double)
   *
   * @param d     Decoded VTK Dataarray plus type tag.
   * @param desc  Fortran descriptor whose base-address will receive the data.
   */
  template<typename T>
  void allocate_and_convert(const DecodedBuffer& d, CFI_cdesc_t* desc);
};

extern "C" {
  /**
   * @brief C-interface constructor for the C++ VTI object.
   *
   * @param vti_path Path to VTI file
   * @return VTI*       VTI object pointer
   */
  gsl::owner<VTI*> VTI__new(const char* vti_path);

  /**
   * @brief Read an integer DataArray into a Fortran pointer descriptor.
   *
   * @param vti   Previously initialized VTI object with allocated vti_tree.
   * @param name  Name of target attribute.
   * @param desc  Pre-allocated descriptor to be filled by \c CFI_allocate.
   */
  void C_VTI_readDatasetInt(VTI* vti, const char* name, CFI_cdesc_t* desc);

  /**
   * @brief Read a floating-point DataArray into a Fortran pointer descriptor.
   *
   * @param vti   Previously initialized VTI object with allocated vti_tree.
   * @param name  Name of target attribute.
   * @param desc  Pre-allocated descriptor to be filled by \c CFI_allocate.
   */
  void C_VTI_readDatasetReal(VTI* vti, const char* name, CFI_cdesc_t* desc);

  /**
   * @brief Extract grid size, physical extent and origin from a VTI file.
   *
   * @param vti            Previously initialized VTI object with allocated vti_tree.
   * @param cells          Number of cells along x/y/z.
   * @param geom_size      Physical side lengths.
   * @param origin         Origin coordinates.
   * @param labels_desc    Optional labels descriptor (may be nullptr).
   */
  void C_VTI_readGeometry(VTI* vti,
                          int* cells,
                          double* geom_size,
                          double* origin,
                          CFI_cdesc_t* labels_desc);

  /**
   * @brief Destroy a VTI instance allocated via VTI__new.
   *
   * @param vti Pointer returned by VTI__new (ignored if nullptr).
   */
  void VTI__delete(gsl::owner<VTI*> vti);
}

#endif
