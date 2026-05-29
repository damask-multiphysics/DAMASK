// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file IO.h
 * @brief C interface to IO.f90 module
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#include <ISO_Fortran_binding.h>

#include <string_view>

extern "C" void F_IO_error(int error_id, CFI_cdesc_t* msg);

class IO {
public:
  inline static void (*fn)(int error_id, CFI_cdesc_t* msg) = &F_IO_error;
  static void error(int error_id, std::string_view msg) {
    CFI_CDESC_T(0) msg_desc_storage;
    auto* msg_desc = reinterpret_cast<CFI_cdesc_t*>(&msg_desc_storage); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    CFI_establish(msg_desc,
                  const_cast<char*>(msg.data()), // NOLINT(cppcoreguidelines-pro-type-const-cast)
                  CFI_attribute_other,
                  CFI_type_char,
                  msg.size(),
                  0,
                  nullptr);
    fn(error_id, msg_desc);
  }
};
