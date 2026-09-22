// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * @file test_VTI.cpp
 * @brief Unit tests for src/grid/VTI.cpp
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#include <gtest/gtest.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <span>
#include <string_view>

#include "../../src/grid/VTI.h"
#include "conftest.h"

struct TempVTIFile {
  TempDirGuard dir;
  std::filesystem::path path;
};

static TempVTIFile write_temp_vti(const std::string& xml) {
  TempVTIFile file{TempDirGuard("damask_vti"), {}};
  file.path = file.dir.path / "test.vti";
  std::ofstream out(file.path);
  if (!out)
    throw std::runtime_error("Failed to create temp VTI file: " + file.path.string());
  out << xml;
  return file;
}

const std::string B64_UC32 = "BAAAAAECAwQAAAAA";             // [1,2,3,4] 32bit uncompressed
const std::string B64_UC64 = "BAAAAAAAAAABAgMEAAAAAAAAAAA="; // [1,2,3,4] 64bit uncompressed

const std::string B64_COMP32 = "AgAAAAQAAAADAAAADAAAAAsAAAB4nGNgZGIGAAAOAAd4nOPi5gEAAEMAIg=="; // [0,1,2,3], [10,11,12] 32bit
                                                                                               // compressed
const std::string B64_COMP64 = "AgAAAAAAAAAEAAAAAAAAAAMAAAAAAAAADAAAAAAAAAALAAAAAAAAAHicY2BkYgYA"
                               "AA4AB3ic4+LmAQAAQwAi"; // [0,1,2,3], [10,11,12] 64bit compressed

const std::vector<uint8_t> EXPECTED_UNCOMPRESSED = {1, 2, 3, 4};
const std::vector<uint8_t> EXPECTED_COMPRESSED = {0, 1, 2, 3, 10, 11, 12};

constexpr std::size_t N_BYTES_PER_WORD_32BIT = 4;
constexpr std::size_t N_BYTES_PER_WORD_64BIT = 8;

TEST(ReadWordTest, ReadsLittleEndian) {
  const std::array<uint8_t, 4> d32 = {0x78, 0x56, 0x34, 0x12};
  const std::array<uint8_t, 8> d64 = {0xF0, 0xDE, 0xBC, 0x9A, 0x78, 0x56, 0x34, 0x12};
  EXPECT_EQ(VTI::read_word(d32), 0x12345678ULL);
  EXPECT_EQ(VTI::read_word(d64), 0x123456789ABCDEF0ULL);
}

TEST(DecodeUncompressedVTI, Uncompressed32Bit) {
  auto out = VTI::decode_uncompressed(B64_UC32, N_BYTES_PER_WORD_32BIT);
  EXPECT_EQ(out, EXPECTED_UNCOMPRESSED);
}

TEST(DecodeUncompressedVTI, Uncompressed64Bit) {
  auto out = VTI::decode_uncompressed(B64_UC64, N_BYTES_PER_WORD_64BIT);
  EXPECT_EQ(out, EXPECTED_UNCOMPRESSED);
}

TEST(DecodeUncompressedVTI, Uncompressed32BitUnderflowError) {
  // specify size 5, only provide 4
  const std::string bad = "BQAAAAECAwQ="; // hex: 05 00 00 00 01 02 03 04
  last_f_io_error_msg().clear();
  EXPECT_THROW(VTI::decode_uncompressed(bad, N_BYTES_PER_WORD_32BIT), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("data block exceeds payload"), std::string::npos);
}

TEST(DecodeCompressedVTI, Compressed32Bit) {
  auto out = VTI::decode_compressed(B64_COMP32, N_BYTES_PER_WORD_32BIT);
  EXPECT_EQ(out, EXPECTED_COMPRESSED);
}

TEST(DecodeCompressedVTI, Compressed64Bit) {
  auto out = VTI::decode_compressed(B64_COMP64, N_BYTES_PER_WORD_64BIT);
  EXPECT_EQ(out, EXPECTED_COMPRESSED);
}

TEST(ReadDatasetRaw, Uncompressed32Bit) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.read_dataset_raw("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, EXPECTED_UNCOMPRESSED);
}

TEST(ReadDatasetRaw, Compressed32Bit) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32" compressor="vtkZLibDataCompressor">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      B64_COMP32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.read_dataset_raw("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, EXPECTED_COMPRESSED);
}

TEST(ReadDatasetRaw, TrailingWhitespaceInBase64) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      B64_UC32 + "   " + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.read_dataset_raw("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, EXPECTED_UNCOMPRESSED);
}

TEST(ReadCellsSizeOrigin, GeometryExtraction) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian">
             <ImageData WholeExtent="0 2 0 3 0 4"
                        Spacing="0.5 1.0 2.0"
                        Origin="10 20 30">
             </ImageData>
           </VTKFile>)";
  std::array<int, 3> cells = {0, 0, 0};
  std::array<double, 3> size = {0, 0, 0};
  std::array<double, 3> org = {0, 0, 0};
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  C_VTI_readGeometry(&vti, cells.data(), size.data(), org.data(), nullptr);
  EXPECT_EQ(cells[0], 2);
  EXPECT_EQ(cells[1], 3);
  EXPECT_EQ(cells[2], 4);
  EXPECT_DOUBLE_EQ(size[0], 1.0); // 0.5 * 2
  EXPECT_DOUBLE_EQ(size[1], 3.0); // 1.0 * 3
  EXPECT_DOUBLE_EQ(size[2], 8.0); // 2.0 * 4
  EXPECT_DOUBLE_EQ(org[0], 10.0);
  EXPECT_DOUBLE_EQ(org[1], 20.0);
  EXPECT_DOUBLE_EQ(org[2], 30.0);
}

TEST(ReadDatasetRaw, ThrowsOnMissingArray) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="other" format="binary">)" +
      B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  last_f_io_error_msg().clear();
  EXPECT_THROW((void)vti.read_dataset_raw("testdata"), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("no DataArray with Name='testdata' found"), std::string::npos);
}

TEST(ParseOnInit, RejectsUnsupportedByteOrder) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="BigEndian">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  last_f_io_error_msg().clear();
  EXPECT_THROW(VTI{file.path.c_str()}, FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("byte_order must be 'LittleEndian'"), std::string::npos);
}

TEST(ParseOnInit, RejectsUnsupportedCompressor) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    compressor="vtkLZ4DataCompressor">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  last_f_io_error_msg().clear();
  EXPECT_THROW(VTI{file.path.c_str()}, FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("compressor is not vtkZLibDataCompressor"), std::string::npos);
}

// Set up parametrization, exposed in tests as TypeParam
// https://google.github.io/googletest/reference/testing.html
template <typename T>
class ReadDatasetInt : public ::testing::Test {};
using ReadDatasetIntTypes = ::testing::Types<int32_t, int64_t>;
TYPED_TEST_SUITE(ReadDatasetInt, ReadDatasetIntTypes);

TYPED_TEST(ReadDatasetInt, ConvertsInt64Input) {
  // CAAAACkAAAAAAAAA -> 41 in Int64
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int64" Name="mydata" format="binary">CAAAACkAAAAAAAAA</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());

  const std::string_view label = "mydata";
  const std::vector<TypeParam> data = vti.read_dataset<TypeParam>(label);
  ASSERT_EQ(data.size(), std::size_t{1});
  EXPECT_EQ(data[0], TypeParam{41});
}
