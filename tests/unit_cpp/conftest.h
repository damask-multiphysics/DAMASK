// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file CLI.h
 * @brief DAMASK C++ unittesting fixtures
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#include <ISO_Fortran_binding.h>

#include <string>
#include <vector>

#include "../../src/IO.h"

// Mock fortran functions (resolve directly so the actual implementations are not pulled from fortran)
namespace {
struct FIOErrorCalled final {};
}

static std::string last_f_io_error_msg_storage;

inline std::string& last_f_io_error_msg() {
  return last_f_io_error_msg_storage;
}

extern "C" {
  static std::vector<std::string> fortran_mock_buffer;

  void F_IO_error(int, CFI_cdesc_t* msg) {
    last_f_io_error_msg() = msg ? std::string(static_cast<char*>(msg->base_addr), msg->elem_len) : "";
    throw FIOErrorCalled{};
  }
  void F_IO_printCppString(CFI_cdesc_t* c_str) {
    fortran_mock_buffer.emplace_back(c_str ? std::string(static_cast<char*>(c_str->base_addr), c_str->elem_len) : "");
  }
  void F_printCompileOptions() {}
  bool IO_redirectedSTDOUT = false;
  bool IO_redirectedSTDERR = false;
}

static bool io_error_called = false;
static int io_error_id = 0;
static std::string io_error_msg;

static void mock_io_error(int error_ID, CFI_cdesc_t* msg) {
  io_error_called = true;
  io_error_id = error_ID;
  std::string msg_str = msg ? std::string(static_cast<char*>(msg->base_addr), msg->elem_len) : "";
  io_error_msg = msg_str;
  last_f_io_error_msg() = msg_str;
  throw FIOErrorCalled{};
}

// Use object lifespan to ensure clean fortran buffer across tests
struct FortranBufferGuard {
  FortranBufferGuard() { fortran_mock_buffer.clear(); last_f_io_error_msg().clear(); }
};

// Use object lifespan to ensure IO mocking across tests
struct IOMockGuard {
  decltype(IO::fn) old_fn;
  IOMockGuard() : old_fn(IO::fn) {
    IO::fn = &mock_io_error;
    io_error_called = false;
    io_error_id = 0;
    io_error_msg.clear();
  }
  ~IOMockGuard() { IO::fn = old_fn; }
};
