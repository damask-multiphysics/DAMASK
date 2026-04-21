#pragma once

#include <ISO_Fortran_binding.h>

#include <string_view>

extern "C" void F_IO_error(int error_ID, CFI_cdesc_t* msg);

class IO {
public:
  inline static void (*fn)(int error_ID, CFI_cdesc_t* msg) = &F_IO_error;
  static void error(int error_ID, std::string_view msg) {
    CFI_CDESC_T(0) msg_desc_storage;
    auto* msg_desc = reinterpret_cast<CFI_cdesc_t*>(&msg_desc_storage);
    CFI_establish(msg_desc, const_cast<char*>(msg.data()), CFI_attribute_other,
                  CFI_type_char, msg.size(), 0, nullptr);
    fn(error_ID, msg_desc);
  }
};
