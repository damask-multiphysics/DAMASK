// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file utils.cpp
 * @brief Shared C++ helper/utility functions for DAMASK solvers
 *
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#ifdef BOOST

#include <vector>
#include <array>
#include <cstddef>
#include <cstring>
#include <string>
#include <cstdint> // IWYU pragma: keep
#include <string_view>
#include <stdexcept>
#include <algorithm>
#include <ISO_Fortran_binding.h>

#include "utils.h"

namespace utils {

template <typename T>
void fortran_array_from_vector(const std::vector<T>& vec, CFI_cdesc_t* array) {
  if (array->base_addr)
    CFI_deallocate(array);
  if (!vec.empty()) {
    const CFI_index_t lower_bound = 1;
    const CFI_index_t extent = static_cast<CFI_index_t>(vec.size());
    if (CFI_allocate(array, &lower_bound, &extent, 0) != CFI_SUCCESS)
      throw std::runtime_error("CFI allocation failed");
    std::memcpy(array->base_addr, vec.data(), vec.size() * sizeof(T));
  }
}
template void fortran_array_from_vector<std::int32_t>(const std::vector<std::int32_t>&, CFI_cdesc_t*);
template void fortran_array_from_vector<std::int64_t>(const std::vector<std::int64_t>&, CFI_cdesc_t*);
template void fortran_array_from_vector<double>(const std::vector<double>&, CFI_cdesc_t*);

void fortran_array_from_vector_str(const std::vector<std::string>& vec,
                                   CFI_cdesc_t* array,
                                   const std::size_t char_len) {
  if (array->base_addr)
    CFI_deallocate(array);
  if (!vec.empty()) {
    const CFI_index_t lower_bound = 1;
    const CFI_index_t extent = static_cast<CFI_index_t>(vec.size());
    if (CFI_allocate(array, &lower_bound, &extent, char_len) != CFI_SUCCESS)
      throw std::runtime_error("CFI allocation failed");
    for (std::size_t i = 0; i < vec.size(); ++i) {
      const std::array<CFI_index_t, 1> sub = {static_cast<CFI_index_t>(i + 1)};
      auto* dst = static_cast<char*>(CFI_address(array, sub.data()));
      std::memset(dst, ' ', char_len);
      const std::size_t copy_len = std::min<std::size_t>(vec.at(i).size(), char_len);
      std::copy_n(vec.at(i).data(), copy_len, dst);
    }
  }
}

std::string_view view_on_fortran_string(const CFI_cdesc_t* string) {
  if (string->type != CFI_type_char)
    throw std::runtime_error("descriptor does not hold a character string");
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  return {reinterpret_cast<const char*>(string->base_addr), string->elem_len};
}

} // namespace utils

#endif
