#include <gtest/gtest.h>
#include <ISO_Fortran_binding.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

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

const std::string kB64UC32 = "BAAAAAECAwQAAAAA";             // [1,2,3,4] 32bit uncompressed
const std::string kB64UC64 = "BAAAAAAAAAABAgMEAAAAAAAAAAA="; // [1,2,3,4] 64it uncompressed

const std::string kB64Comp32 =
    "AgAAAAQAAAADAAAADAAAAAsAAAB4nGNgZGIGAAAOAAd4nOPi5gEAAEMAIg=="; // [0,1,2,3], [10,11,12] 32bit compressed
const std::string kB64Comp64 = "AgAAAAAAAAAEAAAAAAAAAAMAAAAAAAAADAAAAAAAAAALAAAAAAAAAHicY2BkYgYA"
                               "AA4AB3ic4+LmAQAAQwAi"; // [0,1,2,3], [10,11,12] 64bit compressed

const std::vector<uint8_t> kExpectedUncompressed = {1, 2, 3, 4};
const std::vector<uint8_t> kExpectedCompressed = {0, 1, 2, 3, 10, 11, 12};

constexpr std::size_t kNBytesPerWord32bit = 4;
constexpr std::size_t kNBytesPerWord64bit = 8;

TEST(ReadWordTest, ReadsLittleEndian) {
  const std::array<uint8_t, 4> d32 = {0x78, 0x56, 0x34, 0x12};
  const std::array<uint8_t, 8> d64 = {0xF0, 0xDE, 0xBC, 0x9A, 0x78, 0x56, 0x34, 0x12};
  EXPECT_EQ(VTI::read_word(d32.data(), kNBytesPerWord32bit), 0x12345678ULL);
  EXPECT_EQ(VTI::read_word(d64.data(), kNBytesPerWord64bit), 0x123456789ABCDEF0ULL);
}

TEST(DecodeUncompressedVTI, Uncompressed32Bit) {
  auto out = VTI::decode_uncompressed_vti(kB64UC32, kNBytesPerWord32bit);
  EXPECT_EQ(out, kExpectedUncompressed);
}

TEST(DecodeUncompressedVTI, Uncompressed64Bit) {
  auto out = VTI::decode_uncompressed_vti(kB64UC64, kNBytesPerWord64bit);
  EXPECT_EQ(out, kExpectedUncompressed);
}

TEST(DecodeUncompressedVTI, Uncompressed32BitUnderflowError) {
  // specify size 5, only provide 4
  const std::string bad = "BQAAAAECAwQ="; // hex: 05 00 00 00 01 02 03 04
  last_f_io_error_msg().clear();
  EXPECT_THROW(VTI::decode_uncompressed_vti(bad, kNBytesPerWord32bit), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("data block exceeds payload"), std::string::npos);
}

TEST(DecodeCompressedVTI, Compressed32Bit) {
  auto out = VTI::decode_compressed_vti(kB64Comp32, kNBytesPerWord32bit);
  EXPECT_EQ(out, kExpectedCompressed);
}

TEST(DecodeCompressedVTI, Compressed64Bit) {
  auto out = VTI::decode_compressed_vti(kB64Comp64, kNBytesPerWord64bit);
  EXPECT_EQ(out, kExpectedCompressed);
}

TEST(ParseCellDataArray, Uncompressed32Bit) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      kB64UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.parse_cell_data_array("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, kExpectedUncompressed);
}

TEST(ParseCellDataArray, Compressed32Bit) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32" compressor="vtkZLibDataCompressor">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      kB64Comp32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.parse_cell_data_array("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, kExpectedCompressed);
}

TEST(ParseCellDataArray, TrailingWhitespaceInBase64) {
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      kB64UC32 + "   " + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  auto vtk_array = vti.parse_cell_data_array("mydata");
  EXPECT_EQ(vtk_array.vtk_type, "Int32");
  EXPECT_EQ(vtk_array.raw_bytes, kExpectedUncompressed);
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

TEST(ParseCellDataArray, ThrowsOnMissingArray) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="other" format="binary">)" +
      kB64UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  last_f_io_error_msg().clear();
  EXPECT_THROW((void)vti.parse_cell_data_array("testdata"), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("no DataArray with Name='testdata' found"), std::string::npos);
}

TEST(ParseCellDataArray, RejectsUnsupportedByteOrder) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="BigEndian">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      kB64UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  last_f_io_error_msg().clear();
  EXPECT_THROW((void)vti.parse_cell_data_array("mydata"), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("byte_order must be 'LittleEndian'"), std::string::npos);
}

TEST(ParseCellDataArray, RejectsUnsupportedCompressor) {
  IOMockGuard io_mock;
  std::string xml =
      R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    compressor="vtkLZ4DataCompressor">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" +
      kB64UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
  auto file = write_temp_vti(xml);
  VTI vti(file.path.c_str());
  last_f_io_error_msg().clear();
  EXPECT_THROW((void)vti.parse_cell_data_array("mydata"), FIOErrorCalled);
  EXPECT_NE(last_f_io_error_msg().find("compressor is not vtkZLibDataCompressor"), std::string::npos);
}

// Set up parametrization, exposed in tests as TypeParam
// https://google.github.io/googletest/reference/testing.html
template <typename T> class ReadDatasetInt : public ::testing::Test {};
using ReadDatasetIntTypes = ::testing::Types<int32_t, int64_t>;
TYPED_TEST_SUITE(ReadDatasetInt, ReadDatasetIntTypes);

TYPED_TEST(ReadDatasetInt, ConvertsInt64Input) {
  // CAAAACkAAAAAAAAA -> 41 in Int64, + padding equals 42
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

  CFI_CDESC_T(1) desc_storage;
  // https://github.com/gcc-mirror/gcc/blob/master/libgfortran/ISO_Fortran_binding.h#L77
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  CFI_cdesc_t* desc = reinterpret_cast<CFI_cdesc_t*>(&desc_storage);
  std::array<CFI_index_t, 1> extents{};
  CFI_type_t dtype = std::is_same_v<TypeParam, int32_t> ? CFI_type_int32_t : CFI_type_int64_t;
  int rc = CFI_establish(desc, nullptr, CFI_attribute_allocatable, dtype, sizeof(TypeParam), 1, extents.data());
  ASSERT_EQ(rc, CFI_SUCCESS);

  vti.read_dataset_int("mydata", desc);
  std::span<TypeParam> data(static_cast<TypeParam*>(desc->base_addr), static_cast<std::size_t>(desc->dim[0].extent));
  ASSERT_EQ(data.size(), std::size_t{1});
  EXPECT_EQ(data[0], TypeParam{42});
  EXPECT_EQ(CFI_deallocate(desc), CFI_SUCCESS);
}
