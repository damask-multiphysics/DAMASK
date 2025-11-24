#include <gtest/gtest.h>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <system_error>

#ifdef BOOST

#include "../../CLI.h"

namespace fs = std::filesystem;

// Mock fortran functions (resolve directly so the actual implementations are not pulled from fortran)
extern "C" {
  static std::vector<std::string> fortran_mock_buffer;
  void F_IO_printCppString(const char* c_str) { fortran_mock_buffer.emplace_back(c_str ? c_str : ""); }
  void F_printCompileOptions() {}
  bool IO_redirectedSTDOUT = false;
  bool IO_redirectedSTDERR = false;
}

// Use object lifespan to ensure clean fortran buffer across tests
struct FortranBufferGuard {
  FortranBufferGuard()  { fortran_mock_buffer.clear(); }
  ~FortranBufferGuard() { fortran_mock_buffer.clear(); }
};

// Fixture to ensure we go back to the testdir after a test changes our cwd
class CwdGuard : public ::testing::Test {
protected:
  FortranBufferGuard prints_;
  fs::path saved_cwd_;
  void SetUp() override { saved_cwd_ = fs::current_path(); }
  void TearDown() override { std::error_code ec; fs::current_path(saved_cwd_, ec); }
};

TEST(TestSetup, ValidateFortranMockbuffer) {
  FortranBufferGuard guard;
  FortranStream stream;
  stream << "Test string" << std::endl;
  ASSERT_EQ(fortran_mock_buffer.size(), 1u);
  EXPECT_EQ(fortran_mock_buffer.back(), "Test string\n");
}

TEST_F(CwdGuard, SimpleInitialization) {
  int mpi_world_rank = 0;
  const char* argv_literals[] = {
    "dummysolver",
    "--geometry", "=geom.vti", // Test for trailing equal until 4.0
    "--loadcase", "load.yaml",
    "--materialconfig", "material.yaml"
  };
  int argc = static_cast<int>(std::size(argv_literals));
  std::vector<char*> argv;
  argv.reserve(argc);
  for (auto* s : argv_literals)
    argv.push_back(const_cast<char*>(s));

  CLI cli(&argc, argv.data(), &mpi_world_rank);

  EXPECT_EQ(cli.geom_path, fs::path("geom.vti"));
  EXPECT_EQ(cli.loadfile_path, fs::path("load.yaml"));
  EXPECT_EQ(cli.material_path, fs::path("material.yaml"));
  EXPECT_EQ(cli.jobname, "geom_load_material");
  EXPECT_FALSE(cli.uuid.empty());
  EXPECT_FALSE(fortran_mock_buffer.empty());
}

TEST_F(CwdGuard, InitializationWorkingdirNumerics) {
  int mpi_world_rank = 0;

  fs::path test_workdir = fs::temp_directory_path() / "damask_cli_test_wd";
  std::error_code ec;
  fs::create_directories(test_workdir, ec);

  const char* argv_literals[] = {
    "dummysolver",
    "-g", "geom.vti",
    "-l", "load.yaml",
    "-m", "material.yaml",
    "-n", "numerics.yaml",
    "-w", test_workdir.c_str()
  };
  int argc = static_cast<int>(std::size(argv_literals));
  std::vector<char*> argv;
  argv.reserve(argc);
  for (auto* s : argv_literals)
    argv.push_back(const_cast<char*>(s));

  CLI cli(&argc, argv.data(), &mpi_world_rank);

  EXPECT_EQ(cli.numerics_path, fs::path("numerics.yaml"));
  EXPECT_EQ(cli.jobname, "geom_load_material_numerics");
  EXPECT_EQ(fs::current_path(), test_workdir);

  fs::remove_all(test_workdir, ec);
}

TEST_F(CwdGuard, InitializationRestart) {
  int mpi_world_rank = 0;
  const char* argv_literals[] = {
    "dummysolver",
    "--geometry", "geom.vti",
    "--loadcase", "load.yaml",
    "--materialconfig", "material.yaml",
    "--restart", "7"
  };
  int argc = static_cast<int>(std::size(argv_literals));
  std::vector<char*> argv;
  argv.reserve(argc);
  for (auto* s : argv_literals)
    argv.push_back(const_cast<char*>(s));

  CLI cli(&argc, argv.data(), &mpi_world_rank);
  EXPECT_EQ(cli.jobname, "geom_load_material");
  EXPECT_EQ(cli.restart_inc, 7);
}

#endif
