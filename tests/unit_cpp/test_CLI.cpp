// SPDX-License-Identifier: AGPL-3.0-or-later
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <span>
#include <string>
#include <vector>
#include <system_error>
#include <initializer_list>

#include "../../src/CLI.h"
#include "conftest.h"

namespace fs = std::filesystem;

// Fixture to ensure we go back to the testdir after a test changes our cwd
class CwdGuard : public ::testing::Test {
private:
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
  std::vector<const char*> argv = {
    "dummysolver",
    "--geometry", "=geom.vti", // Test for trailing equal until 4.0
    "--loadcase", "load.yaml",
    "--materialconfig", "material.yaml"
  };
  auto args = std::span(argv.data(),std::size(argv));

  CLI cli(args, &mpi_world_rank);

  EXPECT_EQ(cli.geom_path, fs::path("geom.vti"));
  EXPECT_EQ(cli.loadfile_path, fs::path("load.yaml"));
  EXPECT_EQ(cli.material_path, fs::path("material.yaml"));
  EXPECT_EQ(cli.jobname, "geom_load_material");
  EXPECT_EQ(cli.restart_inc, -1);
  EXPECT_FALSE(cli.uuid.empty());
  EXPECT_FALSE(fortran_mock_buffer.empty());
}

TEST_F(CwdGuard, InitializationWorkingdirNumerics) {
  int mpi_world_rank = 0;

  fs::path test_workdir = fs::temp_directory_path() / "damask_cli_test_wd";
  std::error_code ec;
  fs::create_directories(test_workdir, ec);

  std::vector<const char*> argv = {
    "dummysolver",
    "-g", "geom.vti",
    "-l", "load.yaml",
    "-m", "material.yaml",
    "-n", "numerics.yaml",
    "-w", test_workdir.c_str()
  };
  auto args = std::span(argv.data(),std::size(argv));

  CLI cli(args, &mpi_world_rank);

  EXPECT_EQ(cli.numerics_path, fs::path("numerics.yaml"));
  EXPECT_EQ(cli.jobname, "geom_load_material_numerics");
  EXPECT_EQ(fs::current_path(), test_workdir);

  fs::remove_all(test_workdir, ec);
}

TEST_F(CwdGuard, InitializationRestart) {
  int mpi_world_rank = 0;
  std::vector<const char*> argv = {
    "dummysolver",
    "--geometry", "geom.vti",
    "--loadcase", "load.yaml",
    "--materialconfig", "material.yaml",
    "--restart", "7"
  };
  auto args = std::span(argv.data(),std::size(argv));

  CLI cli(args, &mpi_world_rank);
  EXPECT_EQ(cli.jobname, "geom_load_material");
  EXPECT_EQ(cli.restart_inc, 7);
}

TEST_F(CwdGuard, HelpExitsWithoutError) {
  int mpi_world_rank = 0;
  std::vector<const char*> argv = {
    "dummysolver",
    "-h"
  };
  auto args = std::span(argv.data(),std::size(argv));
  EXPECT_EXIT({
      CLI cli(args, &mpi_world_rank);
      (void)cli;
    }, ::testing::ExitedWithCode(0), ""
  );
}

TEST_F(CwdGuard, InitializationRestartNegative) {
  int mpi_world_rank = 0;
  for (const char* restart_val : {"-1", "-2"}) {
    std::vector<const char*> argv = {
      "dummysolver",
      "--geometry", "geom.vti",
      "--loadcase", "load.yaml",
      "--materialconfig", "material.yaml",
      "--restart", restart_val
    };
    auto args = std::span(argv.data(),std::size(argv));

    EXPECT_THROW(CLI cli(args, &mpi_world_rank), FIOErrorCalled);
    EXPECT_EQ(last_f_io_error_msg(), std::string("invalid value for --restart: ") + restart_val);
  }
}

TEST_F(CwdGuard, InitializationMissingLoadcaseValue) {
  int mpi_world_rank = 0;
  IOMockGuard io_mock;

  std::vector<const char*> argv = {
    "dummysolver",
    "-g", "20grains16x16x16.vti",
    "-m", "material.yaml"
  };
  auto args = std::span(argv.data(),std::size(argv));

  EXPECT_THROW(CLI cli(args, &mpi_world_rank), FIOErrorCalled);
  EXPECT_EQ(io_error_id, 610);
  EXPECT_NE(io_error_msg.find("the option '--loadcase' is required but missing"), std::string::npos);
}

TEST_F(CwdGuard, InitializationMissingArgumentAfterFlag) {
  int mpi_world_rank = 0;
  IOMockGuard io_mock;

  std::vector<const char*> argv = {
    "dummysolver",
    "-g", "20grains16x16x16.vti",
    "-m", "material.yaml",
    "-l"
  };
  auto args = std::span(argv.data(),std::size(argv));

  EXPECT_THROW(CLI cli(args, &mpi_world_rank), FIOErrorCalled);
  EXPECT_EQ(io_error_id, 610);
  EXPECT_NE(io_error_msg.find("the required argument for option '--loadcase' is missing"), std::string::npos);
}

TEST_F(CwdGuard, InitializationRestartFromZero) {
  int mpi_world_rank = 0;
  std::vector<const char*> argv = {
    "dummysolver",
    "--geometry", "geom.vti",
    "--loadcase", "load.yaml",
    "--materialconfig", "material.yaml",
    "--restart", "0"
  };
  auto args = std::span(argv.data(),std::size(argv));

  CLI cli(args, &mpi_world_rank);
  EXPECT_EQ(cli.restart_inc, 0);
  bool found = false;
  for (const auto& line : fortran_mock_buffer) {
    if (line.find("Restart increment:  0") != std::string::npos) {
      found = true;
      break;
    }
  }
  EXPECT_TRUE(found);
}
