/**
 * @file CLI.h
 * @brief DAMASK C++ commandline interface for PETSc-based solvers
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#ifndef CLI_H
#define CLI_H

#ifdef BOOST

#include <stdio.h>
#include <string>
#include <filesystem>

namespace boost { namespace program_options { class options_description; } }
namespace po = boost::program_options;
namespace fs = std::filesystem;
using namespace std;

extern "C" {
  /**
   * Print a C null‑terminated string through Fortran's `write` method.
   *
   * @param[in] c_str  Pointer to the first character of C-string (must be 0‑terminated).
   */
  void F_IO_printCppString(const char* c_str);

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
        F_IO_printCppString(buffer.c_str());
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
  int          restart_inc = 0;
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
  CLI(int* argc, char* argv[], int* worldrank);

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
#endif // CLI_H
