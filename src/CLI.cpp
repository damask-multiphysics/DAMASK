// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file CLI.cpp
 * @brief DAMASK C++ commandline interface for PETSc-based solvers
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#ifdef BOOST

#include "CLI.h"
#include <ISO_Fortran_binding.h>
#include <pwd.h>
#include <unistd.h>
#include <array>
#include <boost/asio/ip/host_name.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/system/detail/error_code.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 108600
#include <boost/uuid/basic_random_generator.hpp>
#endif
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <utility>
#include "petscversion.h"

namespace fs = std::filesystem;

std::string IO_color(const std::array<int, 3>& rgb) {
  if (isatty(STDOUT_FILENO) && !IO_redirectedSTDOUT) {
    return "\033[38;2;" + std::to_string(rgb.at(0)) + ";" +
        std::to_string(rgb.at(1)) + ";" +
        std::to_string(rgb.at(2)) + "m";
  }
  return {};
}

std::string IO_color_reset() {
  return isatty(STDOUT_FILENO) && !IO_redirectedSTDOUT ? "\033[0m" : "";
}

CLI::CLI(std::span<const char*> args, int* worldrank) {
  std::string arg_geom, arg_load, arg_material, arg_numerics, arg_wd;
  std::string arg_jobname;
  int arg_rs = -1;

  init_print();

  po::options_description flags(" Valid command line flags");
  flags.add_options()
    ("help,h", "show help")
    ("geometry,g",         po::value<std::string>(&arg_geom)->value_name("[ --geom ]"),         "geometry file")
    ("loadcase,l",         po::value<std::string>(&arg_load)->value_name("[ --load ]"),         "load case")
    ("materialconfig,m",   po::value<std::string>(&arg_material)->value_name("[ --material ]"), "material config")
    ("numericsconfig,n",   po::value<std::string>(&arg_numerics)->value_name("[ --numerics ]"), "numerics config")
    ("jobname,j",          po::value<std::string>(&arg_jobname)->value_name("[ --job ]"),       "job name")
    ("workingdirectory,w", po::value<std::string>(&arg_wd)->value_name("[ --wd ]"),             "working directory")
    ("wd",                 po::value<std::string>(&arg_wd),                                     "alias")
#if defined(GRID) || defined(TEST)
    ("restart,r",          po::value<int>(&arg_rs)->value_name("[ --rs ]"),                     "restart increment")
    ("rs",                 po::value<std::string>(&arg_wd),                                     "alias")
#endif
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(int(args.size()), args.data(), flags), vm);
  po::notify(vm);

  /**
   * @brief Helper method to remove leading equal from argument string.
   *
   * Required for backward compatibility, remove for DAMASK 4.0
   *
   * @param[in] path_str
   */
  auto remove_leading_equal = [](const std::string& arg) -> std::string {
    return arg.length() > 0 && arg.at(0) == '=' ? arg.substr(1):arg;
  };

  /**
   * @brief Get stem for path.
   * @param[in] path_str
   */
  auto stem = [](const std::string& path_str) -> std::string {
    fs::path p(path_str);
    return p.stem().string();
  };

  /* Get username */
  auto get_username = []() -> std::string {
    struct passwd *pw = getpwuid(getuid());
    return std::string(pw != nullptr ? pw->pw_name:"n/a (Error getting username)");
  };

  /**
   * @brief Generate an UUID.
   */
  auto generate_uuid = []() -> std::string {
    return boost::lexical_cast<std::string>(boost::uuids::random_generator()());
  };

  if (vm.count("help") || args.size() == 1) {
    CLI::help_print(flags);
    std::exit(0);
  }

  geom_path = remove_leading_equal(arg_geom);
  loadfile_path = remove_leading_equal(arg_load);
  material_path = remove_leading_equal(arg_material);

  if (!arg_numerics.empty())
    numerics_path = remove_leading_equal(arg_numerics);

  if (!arg_jobname.empty())
    jobname = arg_jobname;
  else {
    jobname  = stem(geom_path) + "_" + stem(loadfile_path) + "_" + stem(material_path);
    if (!arg_numerics.empty())
      jobname += "_" + stem(numerics_path);
  }

  if (!arg_wd.empty()) {
    fs::current_path(remove_leading_equal(arg_wd));
  }

  if (arg_rs != -1)
    restart_inc = arg_rs;

  if (*worldrank==0){
    uuid = generate_uuid();
  }

  boost::system::error_code ec;
  std::string hostname = boost::asio::ip::host_name(ec);
  if (ec) std::runtime_error("Boost hostname collection error: " + ec.message());

  cout << " Host name: " << hostname << std::endl;
  cout << " User name: " << get_username() << std::endl << std::endl;
  cout << " Command line call:  ";
  for (auto& arg : args){
    cout << arg << " ";
  }
  cout << std::endl;
  cout << " Working directory:  " << fs::current_path().string() << std::endl;
  cout << " Geometry:           " << geom_path << std::endl;
  cout << " Load case:          " << loadfile_path << std::endl;
  cout << " Material config:    " << material_path << std::endl;
  if (vm.count("numerics")) {
    cout << " Numerics config:  " << numerics_path << std::endl;
  }
  cout << " Job name:           " << jobname << std::endl;
  cout << " Job ID:             " << uuid << std::endl;
  if (restart_inc > 0) {
    cout << " Restart increment:  " << restart_inc << std::endl;
  }
}

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
void CLI::init_print() {
#ifdef DEBUG
  std::array<int,3> red = {255,0,0};
  cout << IO_color(red);
  cout << "debug version - debug version - debug version - debug version - debug version" << std::endl;
#else
  std::array<int,3> DAMASK_blue = {67,128,208};
  cout << IO_color(DAMASK_blue);
#endif
  cout << R"(
     _/_/_/      _/_/    _/      _/    _/_/      _/_/_/  _/    _/    _/_/_/
    _/    _/  _/    _/  _/_/  _/_/  _/    _/  _/        _/  _/            _/
   _/    _/  _/_/_/_/  _/  _/  _/  _/_/_/_/    _/_/    _/_/          _/_/
  _/    _/  _/    _/  _/      _/  _/    _/        _/  _/  _/            _/
 _/_/_/    _/    _/  _/      _/  _/    _/  _/_/_/    _/    _/    _/_/_/

)";

#ifdef GRID
  std::array<int,3> DAMASK_green = {123,207,68};
  cout << IO_color(DAMASK_green);
  cout << " Grid solver" << std::endl << std::endl;
#elif defined(MESH)
  std::array<int,3> DAMASK_orange = {230,150,68};
  cout << IO_color(DAMASK_orange);
  cout << " Mesh solver" << std::endl << std::endl;
#endif

#ifdef DEBUG
  cout << IO_color(red);
  cout << " debug version - debug version - debug version - debug version - debug version" << std::endl << std::endl;
#endif
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  cout << IO_color_reset();
  cout << " F. Roters et al., Computational Materials Science 158:420–478, 2019" << std::endl
       << " https://doi.org/10.1016/j.commatsci.2018.04.030" << std::endl << std::endl;

#if PETSC_VERSION_MAJOR != 3 || PETSC_VERSION_MINOR < PETSC_MINOR_MIN || PETSC_VERSION_MINOR > PETSC_MINOR_MAX
#error "--  UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION --- UNSUPPORTED PETSc VERSION ---"
#else
  cout << " S. Balay et al., PETSc/TAO User Manual Revision " << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR << std::endl;
#if   PETSC_VERSION_MINOR == 19
cout << " https://doi.org/10.2172/1968587" << endl;
#elif PETSC_VERSION_MINOR == 20
cout << " https://doi.org/10.2172/2205494" << endl;
#elif PETSC_VERSION_MINOR == 21
cout << " https://doi.org/10.2172/2337606" << endl;
#elif PETSC_VERSION_MINOR == 22
cout << " https://doi.org/10.2172/2476320" << endl;
#elif PETSC_VERSION_MINOR == 23
cout << " https://doi.org/10.2172/2565610" << endl;
#elif PETSC_VERSION_MINOR == 24
cout << " https://doi.org/10.2172/2998643" << endl;
#endif
cout << endl;
#endif

  cout << " Version: " << DAMASK_VERSION << std::endl << std::endl;
  cout << " Compiled with:";
#if defined(__clang__)
  cout << " Clang version " << __clang_major__ << "." << __clang_minor__ << "."
       << __clang_patchlevel__ << std::endl;
#elif defined(__GNUC__)
  cout << " GCC version " << __GNUC__ << "." << __GNUC_MINOR__ << "."
       << __GNUC_PATCHLEVEL__ << std::endl;
#elif defined(__INTEL_COMPILER)
  cout << " Intel Compiler version "  << __INTEL_COMPILER << "."
       << __INTEL_COMPILER_UPDATE << std::endl;
#endif
  F_printCompileOptions();
  cout << std::endl;
  cout << " PETSc version: " << PETSC_VERSION_MAJOR << "."
                             << PETSC_VERSION_MINOR << "."
                             << PETSC_VERSION_SUBMINOR << std::endl << std::endl;

  cout << " Compiled at: " << __DATE__ << " at " << __TIME__ << std::endl << std::endl;
}

void CLI::help_print(const po::options_description& flags) {
  cout <<
R"(
 #######################################################################
 DAMASK Command Line Interface:
 Düsseldorf Advanced Material Simulation Kit with PETSc-based solvers
 #######################################################################

)"
  << flags <<
R"(

 -----------------------------------------------------------------------
 Mandatory flags:

  --geom GEOMFILE
)"
#if defined(GRID)
R"(
       Relative or absolute path to a VTK image data file (*.vti)
       with mandatory "material" field variable.
)"
#elif defined(MESH)
R"(
       Relative or absolute path to a Gmsh file (*.msh)
       with definitions of physical groups/tags for material IDs
       and boundary conditions.
)"
#endif
R"(
  --load LOADFILE
       Relative or absolute path to a load case definition
       in YAML format.

  --material MATERIALFILE
       Relative or absolute path to a material configuration
       in YAML format.

 -----------------------------------------------------------------------
 Optional flags:

  --numerics NUMERICSFILE
       Relative or absolute path to a numerics configuration
       in YAML format.

  --jobname JOBNAME
       Job name, defaults to GEOM_LOAD_MATERIAL[_NUMERICS].

  --workingdir WORKINGDIRECTORY
       Working directory, defaults to current directory and
       serves as base directory of relative paths.
)"
#if defined(GRID) || defined(TEST)
R"(
  --restart N
       Restart simulation from given increment.
       Read in increment N and, based on this, continue with
       calculating increments N+1, N+2, ...
       Requires restart information for increment N to be present in
       JOBNAME_restart.hdf5 and will append subsequent results to
       existing file JOBNAME.hdf5.
)"
#endif
R"(
 -----------------------------------------------------------------------
 Help:

  --help
       Display help and exit.

 #######################################################################
)" << std::endl;
}


extern "C" {

  /**
   * @brief C-interface constructor for the C++ CLI object.
   *
   * This is the entry point called from Fortran through `ISO_C_BINDING` to
   * create a new CLI instance that parses the command line.
   *
   * Ownership of the returned pointer is transferred to the Fortran side,
   * which automatically deallocates it when the CLI.f90 module goes out of scope.
   *
   * @param[in,out] argc       Argument count
   * @param[in]     argv       Argument vector
   * @param[in]     worldrank  MPI rank
   * @return CLI*              CLI object pointer
   */
  // NOLINTNEXTLINE(bugprone-reserved-identifier)
  CLI* CLI__new(int* argc, const char* argv[], int* worldrank) {
    auto args = std::span(argv,*argc);
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    return new CLI(args, worldrank);
  }

#ifndef OLD_STYLE_C_TO_FORTRAN_STRING

  /**
   * @brief Transfer parsed CLI arguments to Fortran (GNU >= 12 and Intel C-type descriptors).
   *
   * Each non-empty string member from the `CLI` object is copied into a
   * descriptor of type `CFI_cdesc_t` and allocates the required memory.
   *
   * @param[in]  cli       Pointer to a CLI instance
   * @param[out] geom      Geometry path
   * @param[out] load      Load file path
   * @param[out] material  Material path
   * @param[out] numerics  Numerics path
   * @param[out] uuid      UUID string
   * @param[out] jobname   Job name
   * @param[out] restart   Restart increment (only relevant for grid solver)
   * @param[out] stat      Status flag: 0 on success, 1 on allocation/copy failure
   */
  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  void CLI_getParsedArgs (CLI *cli,
                          CFI_cdesc_t *geom,
                          CFI_cdesc_t *load,
                          CFI_cdesc_t *material,
                          CFI_cdesc_t *numerics,
                          CFI_cdesc_t *uuid,
                          CFI_cdesc_t *jobname,
                          int* restart,
                          int* stat) {
    /**
     * @brief Copy a C++ std::string into a Fortran C descriptor.
     *
     * @param[in]  src   Source C++ string
     * @param[out] dest  Destination Fortran descriptor to reallocate
     * @return           status int 0 on success, 1 on failure
     */
  // NOLINTEND(bugprone-easily-swappable-parameters)
    auto put_string = [](const std::string &src, CFI_cdesc_t *dest) -> int {
      const char *src_c = src.c_str();
      size_t len = std::strlen(src_c);
      if (len != 0) {
        if (CFI_allocate(dest, (CFI_index_t *)0, (CFI_index_t *)0, len) == 0) {
          std::memcpy(dest->base_addr, src_c, len);
          return 0;
        }
        return 1;
      }
      return 1;
    };

    const std::pair<std::string, CFI_cdesc_t*> cli_args[] = {
      {cli->geom_path, geom},
      {cli->loadfile_path, load},
      {cli->material_path, material},
      {cli->numerics_path, numerics},
      {cli->jobname, jobname},
      {cli->uuid, uuid}
    };

    for (const auto &arg : cli_args) {
      if (!arg.first.empty() && put_string(arg.first, arg.second) != 0) {
        *stat = 1;
        return;
      }
    }
#ifdef GRID
    *restart = cli->restart_inc;
#else
    (void)restart;
#endif
    *stat = 0;
  }

#else  // OLD_STYLE_C_TO_FORTRAN_STRING

  /**
   * @brief Copy a C++ std::string to a fixed-size Fortran character buffer.
   *
   * Ensures a trailing null terminator is written (so Fortran may treat it as
   * a C string via `C_CHAR`). If the destination buffer is shorter than the source,
   * it will truncate.
   *
   * @param[in]  src   C++ string to copy.
   * @param[out] dest  Pointer to the Fortran buffer
   */
  static void strcopy_to_fortran(const std::string& src,
                                 char* dest) {
    std::strncpy(dest, src.c_str(), src.length() + 1);
  }

  /**
   * @brief Transfer parsed CLI arguments to Fortran (legacy fixed-length strings).
   *
   * Copies non-empty strings from the `CLI` object into the fortran-provided fixed-size buffers.
   *
   * @param[in]  cli       Pointer to a CLI instanc.
   * @param[out] geom      Geometry path buffer
   * @param[out] load      Load file path buffer
   * @param[out] material  Material path buffer
   * @param[out] numerics  Numerics path buffer
   * @param[out] uuid      UUID buffer
   * @param[out] jobname   Jobname buffer
   * @param[out] restart   Restart increment (only relevant for grid solver)
   * @param[out] stat      Status flag: always returns 0 here, failure will be implicit from buffer error
   */
  void CLI_getParsedArgs(CLI* cli,
                         char geom[],
                         char load[],
                         char material[],
                         char numerics[],
                         char uuid[],
                         char jobname[],
                         int* restart,
                         int* stat) {
    if (!cli->geom_path.empty())
      strcopy_to_fortran(cli->geom_path, geom);
    if (!cli->loadfile_path.empty())
      strcopy_to_fortran(cli->loadfile_path, load);
    if (!cli->material_path.empty())
      strcopy_to_fortran(cli->material_path, material);
    if (!cli->numerics_path.empty())
      strcopy_to_fortran(cli->numerics_path, numerics);
    if (!cli->uuid.empty())
      strcopy_to_fortran(cli->uuid, uuid);
    if (!cli->jobname.empty())
      strcopy_to_fortran(cli->jobname, jobname);
#ifdef GRID
    *restart = cli->restart_inc;
#else
    (void)restart;
#endif
    *stat = 0;
  }
#endif
}

#endif
