// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file utils.h
 * @brief Shared C++ helper/utility functions
 *
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#ifdef BOOST

#include <ISO_Fortran_binding.h>
#include <cstddef>
#include <string>
#include <cstdint> // IWYU pragma: keep
#include <string_view>
#include <vector>

namespace utils {

constexpr std::size_t STRLEN = 256;

/**
 * @brief Create a Fortran array and fill it with the content of a vector.
 *
 * @param vec[in] C++ vector.
 * @param array Fortran array.
 */
template <typename T>
void fortran_array_from_vector(const std::vector<T>& vec, CFI_cdesc_t* array);
// Explicit instantiation declarations: the definitions live in utils.cpp.
extern template void fortran_array_from_vector<std::int32_t>(const std::vector<std::int32_t>&, CFI_cdesc_t*);
extern template void fortran_array_from_vector<std::int64_t>(const std::vector<std::int64_t>&, CFI_cdesc_t*);
extern template void fortran_array_from_vector<double>(const std::vector<double>&, CFI_cdesc_t*);

/**
 * @brief Create a Fortran array and fill it with the content of a string vector.
 *
 * @param vec[in] C++ vector of strings.
 * @param array Fortran array of strings.
 */
void fortran_array_from_vector_str(const std::vector<std::string>& vec,
                                   CFI_cdesc_t* array,
                                   const std::size_t char_len = STRLEN);

/**
 * @brief Provide a view on a Fortran string.
 *
 * @param[in] string Fortran string descriptor.
 * @return View on the Fortran string.
 */
std::string_view view_on_fortran_string(const CFI_cdesc_t* string);

} // namespace utils

#endif
