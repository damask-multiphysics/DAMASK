// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file CLI.h
 * @brief DAMASK C++ commandline interface for PETSc-based solvers
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#pragma once

#ifdef BOOST

#include <stdio.h>
#include <ISO_Fortran_binding.h>
#include <ostream>
#include <span>
#include <string>

namespace boost { namespace program_options { class options_description; } }
namespace po = boost::program_options;
using namespace std;

extern "C" {
  /**
   * Print a descriptor-backed string through Fortran's `write` method.
   *
   * @param[in] c_str  Descriptor for the string to print.
   */
  void F_IO_printCppString(CFI_cdesc_t* c_str);

  /** Print Fortran `compiler_options()` string and cmake info. */
  void F_printCompileOptions();

  extern bool IO_redirectedSTDOUT;
  extern bool IO_redirectedSTDERR;
}

/**
 * @class FortranStream
 * @brief Inherited class from `std::ostream`
 *
 * Overrides `overflow` and `sync` to flush stream to Fortran `write` statement.
 * https://en.cppreference.com/w/cpp/io/basic_streambuf.html
 *
 */
class FortranStream : public std::ostream {
  class FortranBuffer : public std::streambuf {
    std::string buffer;
    int overflow(int c) override {
      if (c != EOF) buffer += static_cast<char>(c);
      return c;
    }
    int sync() override {
      if (!buffer.empty()) {
        // use capitalized cdesc_t as macro to create cfi scalar on stack, for storage of string data for Fortran
        // https://www.ibm.com/docs/en/xl-fortran-linux/16.1.1?topic=29113-example-allocatable-pointer-arguments
        CFI_CDESC_T(0) buffer_desc_raw;
        auto* buffer_desc = reinterpret_cast<CFI_cdesc_t*>(&buffer_desc_raw);
        CFI_establish(buffer_desc, buffer.data(), CFI_attribute_other, CFI_type_char, buffer.size(), 0, nullptr);
        F_IO_printCppString(buffer_desc);
        buffer.clear();
      }
      return 0;
    }
  };
  FortranBuffer buf_;
public:
  FortranStream() : std::ostream(&buf_) {}
};

/**
 * Parse command‑line arguments.
 *
 * Instance is created and managed by the Fortran side.
 * Class attributes are passed to Fortran via external CLI_getParsedArgs method.
 */
class CLI {
public:
  int          restart_inc = -1;
  std::string  geom_path, loadfile_path,
               material_path, numerics_path,
               jobname, uuid;

  FortranStream cout{};

  /**
   * @brief Construct CLI instance.
   * @param[in,out] argc        Argument count
   * @param[in]     argv        Argument vector
   * @param[in]     worldrank   MPI rank
   */
  CLI(std::span<const char*> args, int* worldrank);

  /** Print initialization text. */
  void init_print();

  /**
   * @brief Print the help text.
   * @param[in] flags            Boost flags
   * @param[in] include_restart  Print restart help (only grid)
   */
  void help_print(const po::options_description& flags);
};

#endif
