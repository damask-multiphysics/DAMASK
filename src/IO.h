#pragma once

#include <string>
#include <string_view>

extern "C" void F_IO_error(int error_ID, const char* msg);

class IO {
public:
  void (*fn)(int error_ID, const char* msg) = &F_IO_error;
  void error(int error_ID, std::string_view msg) const {
    const std::string tmp(msg);
    fn(error_ID, tmp.c_str());
  }
};
