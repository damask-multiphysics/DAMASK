#include <gtest/gtest.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/version.hpp>
#if BOOST_VERSION >= 108800
#include <boost/uuid/basic_random_generator.hpp>
#endif
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "../../src/grid/VTI.h"
#include "conftest.h"

static std::string write_temp_vti(const std::string& xml) {
    boost::uuids::uuid uuid = boost::uuids::random_generator()();
    std::string filename = boost::uuids::to_string(uuid) + ".vti";
    std::filesystem::path path = std::filesystem::temp_directory_path() / filename;
    std::ofstream out(path);
    if (!out)
      throw std::runtime_error("Failed to create temp VTI file: " + path.string());
    out << xml;
    return path.string();
}

const std::string B64_UC32 = "BAAAAAECAwQAAAAA";                 // [1,2,3,4] 32bit uncompressed
const std::string B64_UC64 = "BAAAAAAAAAABAgMEAAAAAAAAAAA=";     // [1,2,3,4] 64it uncompressed

const std::string B64_COMP32 =
    "AgAAAAQAAAADAAAADAAAAAsAAAB4nGNgZGIGAAAOAAd4nOPi5gEAAEMAIg==";     // [0,1,2,3], [10,11,12] 32bit compressed
const std::string B64_COMP64 =
    "AgAAAAAAAAAEAAAAAAAAAAMAAAAAAAAADAAAAAAAAAALAAAAAAAAAHicY2BkYgYA"
    "AA4AB3ic4+LmAQAAQwAi";                                             // [0,1,2,3], [10,11,12] 64bit compressed

const std::vector<uint8_t> k_expected_uncompressed = {1,2,3,4};
const std::vector<uint8_t> k_expected_compressed = {0,1,2,3,10,11,12};

constexpr std::size_t n_bytes_per_word_32bit = 4;
constexpr std::size_t n_bytes_per_word_64bit = 8;

TEST(ReadWordTest, ReadsLittleEndian) {
    const std::array<uint8_t, 4> d32 = {0x78,0x56,0x34,0x12};
    const std::array<uint8_t, 8> d64 = {0xF0,0xDE,0xBC,0x9A,0x78,0x56,0x34,0x12};
    EXPECT_EQ(VTI::read_word(d32.data(), n_bytes_per_word_32bit), 0x12345678ULL);
    EXPECT_EQ(VTI::read_word(d64.data(), n_bytes_per_word_64bit), 0x123456789ABCDEF0ULL);
}

TEST(DecodeUncompressedVTI, Uncompressed32Bit) {
    auto out = VTI::decode_uncompressed_VTI(B64_UC32, n_bytes_per_word_32bit);
    EXPECT_EQ(out, k_expected_uncompressed);
}

TEST(DecodeUncompressedVTI, Uncompressed64Bit) {
    auto out = VTI::decode_uncompressed_VTI(B64_UC64, n_bytes_per_word_64bit);
    EXPECT_EQ(out, k_expected_uncompressed);
}

TEST(DecodeUncompressedVTI, Uncompressed32BitUnderflowError) {
    // specify size 5, only provide 4
    const std::string bad = "BQAAAAECAwQ="; // hex: 05 00 00 00 01 02 03 04
    last_f_io_error_msg().clear();
    EXPECT_THROW(VTI::decode_uncompressed_VTI(bad, n_bytes_per_word_32bit), FIOErrorCalled);
    EXPECT_NE(last_f_io_error_msg().find("data block exceeds payload"), std::string::npos);
}

TEST(DecodeCompressedVTI, Compressed32Bit) {
    auto out = VTI::decode_compressed_VTI(B64_COMP32, n_bytes_per_word_32bit);
    EXPECT_EQ(out, k_expected_compressed);
}

TEST(DecodeCompressedVTI, Compressed64Bit) {
    auto out = VTI::decode_compressed_VTI(B64_COMP64, n_bytes_per_word_64bit);
    EXPECT_EQ(out, k_expected_compressed);
}

TEST(ParseCellDataArray, Uncompressed32Bit) {
    std::string xml =
        R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" + B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
    auto file = write_temp_vti(xml);
    VTI vti(file.c_str());
    auto vtk_array = vti.parse_cell_data_array("mydata");
    EXPECT_EQ(vtk_array.vtk_type, "Int32");
    EXPECT_EQ(vtk_array.raw_bytes, k_expected_uncompressed);
}

TEST(ParseCellDataArray, Compressed32Bit) {
    std::string xml =
        R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32" compressor="vtkZLibDataCompressor">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)" + B64_COMP32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
    auto file = write_temp_vti(xml);
    VTI vti(file.c_str());
    auto vtk_array = vti.parse_cell_data_array("mydata");
    EXPECT_EQ(vtk_array.vtk_type, "Int32");
    EXPECT_EQ(vtk_array.raw_bytes, k_expected_compressed);
}

TEST(ParseCellDataArray, TrailingWhitespaceInBase64) {
    std::string xml =
        R"(<?xml version="1.0"?>
           <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian"
                    header_type="UInt32">
             <ImageData WholeExtent="0 1 0 1 0 1" Spacing="1 1 1" Origin="0 0 0">
               <Piece Extent="0 1 0 1 0 1">
                 <CellData>
                   <DataArray type="Int32" Name="mydata" format="binary">)"
        + B64_UC32 + "   " + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
    auto file = write_temp_vti(xml);
    VTI vti(file.c_str());
    auto vtk_array = vti.parse_cell_data_array("mydata");
    EXPECT_EQ(vtk_array.vtk_type, "Int32");
    EXPECT_EQ(vtk_array.raw_bytes, k_expected_uncompressed);
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
    std::array<int, 3> cells = {0,0,0};
    std::array<double, 3> size = {0,0,0};
    std::array<double, 3> org  = {0,0,0};
    auto file = write_temp_vti(xml);
    VTI vti(file.c_str());
    C_VTI_readGeometry(&vti, cells.data(), size.data(), org.data(), nullptr);
    EXPECT_EQ(cells[0], 2);
    EXPECT_EQ(cells[1], 3);
    EXPECT_EQ(cells[2], 4);
    EXPECT_DOUBLE_EQ(size[0], 1.0);  // 0.5 * 2
    EXPECT_DOUBLE_EQ(size[1], 3.0);  // 1.0 * 3
    EXPECT_DOUBLE_EQ(size[2], 8.0);  // 2.0 * 4
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
                   <DataArray type="Int32" Name="other" format="binary">)"
        + B64_UC32 + R"(</DataArray>
                 </CellData>
               </Piece>
             </ImageData>
           </VTKFile>)";
    auto file = write_temp_vti(xml);
    VTI vti(file.c_str());
    last_f_io_error_msg().clear();
    EXPECT_THROW((void)vti.parse_cell_data_array("testdata"), FIOErrorCalled);
    EXPECT_NE(last_f_io_error_msg().find("no DataArray with Name='testdata' found"), std::string::npos);
}
