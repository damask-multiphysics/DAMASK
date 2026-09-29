/**
 * @file VTI.h
 * @brief DAMASK VTI reader with beast-based XML parser
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#ifdef BOOST

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "ISO_Fortran_binding.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

namespace pt = boost::property_tree;
namespace fs = std::filesystem;

struct DecodedBuffer {
  std::string vtk_type;
  std::vector<std::uint8_t> raw_bytes;
};

/**
 * @brief Holds the parsed VTI property tree.
 */
class VTI {
public:
  static constexpr std::size_t DIM = 3;

  /**
   * @brief Construct by reading and parsing a VTI file from disk.
   *
   * @param file_path Path to a .vti file.
   * @throws std::runtime_error On I/O or parse/validation errors.
   */
  VTI(const fs::path& file_path);

  /**
   * @brief Read a little-endian word from @p word.
   *
   * @param word  Span holding a single word of size 4 or 8 bytes.
   * @return      Zero-extended 64-bit unsigned integer result.
   */
  static std::uint64_t read_word(std::span<const std::uint8_t> word);

  /**
   * @brief Inflate and concatenate all compressed blocks in a VTI DataArray.
   *
   * @param b64_string  Raw Base-64 string
   * @param n_bytes_per_word Word size (4 or 8).
   * @return            Vector with decoded uncompressed bytes
   */
  static std::vector<std::uint8_t> decode_compressed(std::string_view b64_string, std::size_t n_bytes_per_word);

  /**
   * @brief Decode an uncompressed VTI Dataarray.
   *
   * @param b64_string  Raw Base-64 string.
   * @param n_bytes_per_word Word size (4 or 8).
   * @return            Vector with decoded bytes.
   */
  static std::vector<std::uint8_t> decode_uncompressed(std::string_view b64_string, std::size_t n_bytes_per_word);

  /**
   * @brief Read a DataArray and return the values as a vector.
   *
   * @tparam T    Target type (std::int32_t, std::int64_t, or double).
   * @param label Name of target attribute.
   * @return      Converted values.
   */
  template <typename T>
    requires(std::same_as<T, std::int32_t> || std::same_as<T, std::int64_t> || std::same_as<T, double>)
  std::vector<T> read_dataset(const std::string_view label) const;

  /**
   * @brief Extract grid size, physical extent, origin and cell data labels from a VTI file.
   *
   * @param cells       Number of cells along x/y/z (output).
   * @param geom_size   Physical side lengths (output).
   * @param origin      Origin coordinates (output).
   * @param labels      Cell data labels (output).
   */
  void read_geometry(std::span<int, VTI::DIM> cells,
                     std::span<double, VTI::DIM> geom_size,
                     // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                     std::span<double, VTI::DIM> origin,
                     std::vector<std::string>& labels) const;

  /**
   * @brief Locate a VTK DataArray inside the class VTKFile buffer and return its bytes.
   *
   * @param label Name of the target DataArray.
   * @return      Struct with vtk datatype and the raw decoded bytes.
   */
  DecodedBuffer read_dataset_raw(const std::string_view label) const;

private:
  pt::ptree tree;

  /**
   * @brief Decodes a Base-64 string using Boost.Beast.
   *
   * @param b64 ASCII string with valid Base-64 characters
   * @return Vector with decoded bytes
   */
  static std::vector<std::uint8_t> decode_b64(std::string_view b64);
};

// Explicit instantiation declarations: the definitions live in VTI.cpp.
extern template std::vector<std::int32_t> VTI::read_dataset<std::int32_t>(std::string_view) const;
extern template std::vector<std::int64_t> VTI::read_dataset<std::int64_t>(std::string_view) const;
extern template std::vector<double> VTI::read_dataset<double>(std::string_view) const;

extern "C" {
/**
 * @brief C-interface constructor for the C++ VTI object.
 *
 * @param vti_path Path to VTI file
 * @return VTI*    Owning pointer to the created object; must be released with C_VTI_delete.
 */
VTI* C_VTI_new(const CFI_cdesc_t* vti_path);

/**
 * @brief Read an integer DataArray into a Fortran pointer descriptor.
 *
 * @param vti       Previously initialized VTI object with allocated tree.
 * @param label     Name of target attribute.
 * @param array_out Pre-allocated descriptor to be filled by \c CFI_allocate.
 */
void C_VTI_readDatasetInt(const VTI* vti, const CFI_cdesc_t* label, CFI_cdesc_t* array_out);

/**
 * @brief Read a floating-point DataArray into a Fortran pointer descriptor.
 *
 * @param vti       Previously initialized VTI object with allocated tree.
 * @param label     Name of target attribute.
 * @param array_out Pre-allocated descriptor to be filled by \c CFI_allocate.
 */
void C_VTI_readDatasetReal(const VTI* vti, const CFI_cdesc_t* label, CFI_cdesc_t* array_out);

/**
 * @brief Extract grid size, physical extent and origin from a VTI file.
 *
 * @param vti       Previously initialized VTI object with allocated tree.
 * @param cells     Number of cells along x/y/z.
 * @param geom_size Physical side lengths.
 * @param origin    Origin coordinates.
 * @param array_out Optional labels descriptor (may be nullptr).
 */
void C_VTI_readGeometry(const VTI* vti, int* cells, double* geom_size, double* origin, CFI_cdesc_t* array_out);

/**
 * @brief Destroy a VTI instance allocated via VTI__new.
 *
 * @param vti Owning pointer returned by C_VTI_new (ignored if nullptr).
 */
void C_VTI_delete(VTI* vti);
}

#endif
