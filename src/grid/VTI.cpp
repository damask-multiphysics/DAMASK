/**
 * @file VTI.cpp
 * @brief DAMASK VTI reader with beast-based XML parser
 *
 * @author Daniel Otto de Mentock, Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 * @author Martin Diehl, KU Leuven
 * @copyright
 *   Max‑Planck‑Institut für Nachhaltige Materialien GmbH
 */

#ifdef BOOST

#include <zconf.h>
#include <zlib.h>
#include <algorithm>
#include <array>
#include <ISO_Fortran_binding.h>
#include <cctype>
#include <cstring>
#include <fstream>
#include <functional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <boost/beast/core/detail/base64.hpp>
#include <boost/optional/optional.hpp>
#include <boost/optional/detail/optional_reference_spec.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 108800
#include <boost/beast/core/detail/base64.ipp>
#include <boost/core/addressof.hpp>
#include <boost/iterator/iterator_facade.hpp>
#endif

#include "VTI.h"
#include "../IO.h"

constexpr int kVTKError = 844;

VTI::VTI(const char* file_path) {
  if (!file_path)
    IO::error(kVTKError, "no valid geometry file path supplied");
  this->file_path = file_path;
  std::ifstream f(file_path, std::ios::binary);
  if (!f)
    IO::error(kVTKError, std::string("cannot open file '") + file_path + '\'');
  pt::read_xml(f, vti_tree);
}

uint64_t VTI::read_word(const uint8_t* p, const std::size_t n_bytes_per_word) {
    uint64_t w = 0;
    std::memcpy(&w, p, n_bytes_per_word);
    return w;
};

std::vector<uint8_t> VTI::decode_compressed_VTI(const std::string& b64_string, const std::size_t n_bytes_per_word) {
  const std::vector<uint8_t> decoded_vec = VTI::decode_b64(b64_string);
  const std::span<const uint8_t> decoded(decoded_vec);

  const std::size_t n_header_bytes = 3 * n_bytes_per_word; // VTK generates headers of size 3
  if (decoded.size() < n_header_bytes)
    IO::error(kVTKError, "header for compressed VTI too short");
  std::span<const uint8_t> header = decoded.first(n_header_bytes);
  uint64_t n_blocks               = VTI::read_word(header.subspan(0 * n_bytes_per_word).data(), n_bytes_per_word);
  uint64_t block_uncompressed     = VTI::read_word(header.subspan(1 * n_bytes_per_word).data(), n_bytes_per_word);
  uint64_t last_block_size        = VTI::read_word(header.subspan(2 * n_bytes_per_word).data(), n_bytes_per_word);

  const std::size_t header_size  = n_header_bytes + static_cast<std::size_t>(n_blocks) * n_bytes_per_word;
  if (decoded.size() < header_size)
    IO::error(kVTKError, "missing header for compressed VTI");
  std::span<const uint8_t> c_table = decoded.subspan(n_header_bytes, header_size - n_header_bytes);
  std::vector<std::size_t> block_sizes(n_blocks);
  for (std::size_t i = 0; i < n_blocks; ++i) {
    uint64_t word = VTI::read_word(c_table.subspan(i * n_bytes_per_word).data(), n_bytes_per_word);
    block_sizes[i] = static_cast<std::size_t>(word);
  }

  std::span<const uint8_t> deflated = decoded.subspan(header_size);
  const std::size_t total_uncompressed =
      (n_blocks > 1 ? (n_blocks - 1) * static_cast<std::size_t>(block_uncompressed) : 0) +
      (last_block_size ? static_cast<std::size_t>(last_block_size)
                       : static_cast<std::size_t>(block_uncompressed));
  std::vector<uint8_t> res_vec(total_uncompressed);
  std::span<uint8_t> res(res_vec);
  std::size_t src_offset = 0;
  std::size_t dst_offset = 0;
  for (std::size_t block_idx = 0; block_idx < n_blocks; ++block_idx) {
    if (src_offset + block_sizes[block_idx] > deflated.size())
      IO::error(kVTKError, "invalid VTI, defined size overflowing for block " + std::to_string(block_idx));

    const std::size_t uncompressed_len =
        (block_idx + 1 == n_blocks && last_block_size > 0)
            ? static_cast<std::size_t>(last_block_size)
            : static_cast<std::size_t>(block_uncompressed);
    uLongf dest = static_cast<uLongf>(uncompressed_len);
    int z = uncompress(res.subspan(dst_offset).data(),
                       &dest,
                       deflated.subspan(src_offset).data(),
                       static_cast<uLongf>(block_sizes[block_idx]));
    if (z != Z_OK || dest != static_cast<uLongf>(uncompressed_len))
      IO::error(kVTKError, "zlib inflate failed on block " + std::to_string(block_idx));
    src_offset += block_sizes[block_idx];
    dst_offset += uncompressed_len;
  }
  return res_vec;
}

std::vector<uint8_t> VTI::decode_uncompressed_VTI(const std::string& b64_string, const std::size_t n_bytes_per_word) {
  const std::vector<uint8_t> decoded_vec = VTI::decode_b64(b64_string);
  const std::span<const uint8_t> decoded(decoded_vec);

  std::vector<uint8_t> res_vec;
  std::size_t offset = 0;

  while (offset + n_bytes_per_word <= decoded.size()) {
    const uint64_t n_bytes = VTI::read_word(decoded.subspan(offset).data(), n_bytes_per_word);
    offset += n_bytes_per_word;

    if (n_bytes == 0) break;
    if (offset + n_bytes > decoded.size())
      IO::error(kVTKError, "VTI uncompressed: data block exceeds payload");
    std::span<const uint8_t> payload = decoded.subspan(offset, static_cast<std::size_t>(n_bytes));
    res_vec.insert(res_vec.end(), payload.begin(), payload.end());
    offset += static_cast<std::size_t>(n_bytes);
  }
  return res_vec;
}

void VTI::read_dataset_int(const char* name, CFI_cdesc_t* desc) {
  DecodedBuffer d = parse_cell_data_array(name);
  allocate_and_convert<int>(d, desc);
  if (desc->elem_len == sizeof(int32_t)) {
    std::span<int32_t> data(static_cast<int32_t*>(desc->base_addr),
                            static_cast<std::size_t>(desc->dim[0].extent));
    for (int32_t& v : data) v += 1;
  } else if (desc->elem_len == sizeof(int64_t)) {
    std::span<int64_t> data(static_cast<int64_t*>(desc->base_addr),
                            static_cast<std::size_t>(desc->dim[0].extent));
    for (int64_t& v : data) v += 1;
  }
}

void VTI::read_dataset_real(const char* name, CFI_cdesc_t* desc) {
  DecodedBuffer d = parse_cell_data_array(name);
  allocate_and_convert<double>(d, desc);
}

void VTI::read_geometry(int* cells_ptr,
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        double* geom_size_ptr,
                        double* origin_ptr,
                        CFI_cdesc_t* labels_desc) {
  constexpr std::size_t vec_size = 3;

  auto parse_ints = [](const std::string& s) {
    std::istringstream is(s);
    std::vector<long long> out;
    long long v = 0;
    while (is >> v) out.push_back(v);
    return out;
  };

  auto parse_3_doubles = [&](const std::string& field_name, const std::string& s) {
    std::istringstream is(s);
    std::array<double, 3> out{};
    is >> out[0] >> out[1] >> out[2];
    if (!is)
      IO::error(kVTKError, "bad numeric field for '" + field_name + "' (got '" + s + "')");
    return out;
  };

  auto root = vti_tree.get_child_optional("VTKFile");
  if (!root)
    IO::error(kVTKError, "missing <VTKFile> element");
  const std::string vtkfile_type = get_attr(*root, "type");
  if (vtkfile_type != "ImageData")
    IO::error(kVTKError, "not an ImageData VTK file (type='" + vtkfile_type + "')");
  auto img = root->get_child_optional("ImageData");
  if (!img)
    IO::error(kVTKError, "missing ImageData element");
  std::string dir = get_attr(*img, "Direction");
  if (!dir.empty() && dir != "1 0 0 0 1 0 0 0 1")
    IO::error(kVTKError, "unsupported 'Direction' (got '" + dir + "')");
  const std::string extent_str = get_attr(*img, "WholeExtent");
  if (extent_str.empty())
    IO::error(kVTKError, "missing 'WholeExtent'");
  const std::vector<long long> extent = parse_ints(extent_str);
  if (extent.size() != 2 * vec_size || (extent[0] != 0 || extent[2] != 0 || extent[4] != 0))
    IO::error(kVTKError, "invalid 'WholeExtent' (got '" + extent_str + "')");

  const std::array<double, vec_size> spacing = parse_3_doubles("ImageData@Spacing", get_attr(*img, "Spacing"));
  const std::array<double, vec_size> origin_vec = parse_3_doubles("ImageData@Origin", get_attr(*img, "Origin"));

  std::span<int, vec_size> cells(cells_ptr, vec_size);
  std::span<double, vec_size> geom_size(geom_size_ptr, vec_size);
  std::span<double, vec_size> origin(origin_ptr, vec_size);

  /* modern form for component-wise assignment below (does not work on macOS)
  auto cells_vec = extent | std::views::drop(1) | std::views::stride(2);
  std::copy(cells_vec.begin(), cells_vec.end(), cells.begin());
  */
  cells[0] = int(extent[1]);
  cells[1] = int(extent[3]);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
  cells[2] = int(extent[5]);
  geom_size[0] = spacing[0] * cells[0];
  geom_size[1] = spacing[1] * cells[1];
  geom_size[2] = spacing[2] * cells[2];
  std::copy(origin_vec.begin(), origin_vec.end(), origin.begin());

  if (geom_size[0] <= 0 || geom_size[1] <= 0 || geom_size[2] <= 0)
    IO::error(kVTKError, "one or more entries <= 0 for 'size'");
  if (cells[0] < 1 || cells[1] < 1 || cells[2] < 1)
    IO::error(kVTKError, "one or more entries < 1 for 'cells'");

  if (labels_desc) {
    if (auto cell_data = img->get_child_optional("Piece.CellData")) {
      std::vector<std::string> labels;
      labels.reserve(cell_data->size());
      for (auto& child : *cell_data) {
        if (child.first != "DataArray") continue;
        std::string name = get_attr(child.second, "Name");
        if (name.empty()) continue;
        if (std::find(labels.begin(), labels.end(), name) != labels.end()) {
          IO::error(kVTKError, "repeated label '" + name + '\'');
        }
        labels.push_back(std::move(name));
      }
      if (labels_desc->attribute != CFI_attribute_allocatable &&
          labels_desc->attribute != CFI_attribute_pointer)
        throw std::runtime_error("labels descriptor must be allocatable or pointer");
      if (labels_desc->base_addr)
        CFI_deallocate(labels_desc);
      if (!labels.empty()) {
        constexpr std::size_t char_len = 256;
        const CFI_index_t lb  = 1;
        const CFI_index_t ext = static_cast<CFI_index_t>(labels.size());
        if (CFI_allocate(labels_desc, &lb, &ext, char_len) != CFI_SUCCESS)
          throw std::runtime_error("CFI allocation failed for labels");

        for (std::size_t i = 0; i < labels.size(); ++i) {
          std::array<CFI_index_t, 1> sub = { static_cast<CFI_index_t>(i + 1) };
          char* dst = static_cast<char*>(CFI_address(labels_desc, sub.data()));
          std::memset(dst, ' ', char_len);
          const std::size_t copy_len = std::min<std::size_t>(labels[i].size(), char_len);
          std::copy_n(labels[i].data(), copy_len, dst);
        }
      }
    }
  }
}

DecodedBuffer VTI::parse_cell_data_array(const char* array_name) {
  boost::optional<pt::ptree&> root = vti_tree.get_child_optional("VTKFile");
  if (!root)
    IO::error(kVTKError, "missing <VTKFile> element");
  const std::string type = get_attr(*root, "type");
  if (type != "ImageData")
    IO::error(kVTKError, "VTKFile@type is not ImageData (got '" + type + "')");
  boost::optional<pt::ptree&> image_data = root->get_child_optional("ImageData");
  if (!image_data)
    IO::error(kVTKError, "missing <ImageData> element");
  boost::optional<pt::ptree&> cell_data = image_data->get_child_optional("Piece.CellData");
  if (!cell_data)
    IO::error(kVTKError, "missing <CellData> element");
  boost::optional<pt::ptree&> data_array_node;
  for (auto& child : *cell_data) {
    if (child.first == "DataArray") {
      auto name_attr = get_attr(child.second, "Name");
      if (name_attr == array_name) {
        data_array_node = child.second;
        break;
      }
    }
  }
  if (!data_array_node)
    IO::error(kVTKError, "no DataArray with Name='" + std::string(array_name) + "' found");
  if (get_attr(*data_array_node, "format") != "binary")
    IO::error(kVTKError, "DataArray '" + std::string(array_name) + "' is not binary");
  const std::string vtk_type = get_attr(*data_array_node, "type");
  if (vtk_type.empty())
    IO::error(kVTKError, "DataArray missing 'type' attribute");
  bool compressed = (get_attr(*root, "compressor") == "vtkZLibDataCompressor");
  std::string header_type = get_attr(*root, "header_type");
  if (header_type.empty())
    header_type = "UInt32";
  const std::size_t n_bytes_per_word = (header_type == "UInt64") ? 8 : 4;

  std::string text_content = data_array_node->get_value<std::string>();
  std::vector<uint8_t> raw_bytes = compressed
      ? decode_compressed_VTI(text_content.c_str(), n_bytes_per_word)
      : decode_uncompressed_VTI(text_content.c_str(), n_bytes_per_word);

  return DecodedBuffer{std::move(vtk_type), std::move(raw_bytes)};
}

std::string VTI::get_attr(const pt::ptree& n, const char* key) const {
  if (auto a = n.get_child_optional("<xmlattr>")) {
    if (auto v = a->get_optional<std::string>(key))
      return *v;
  }
  return {};
}

std::vector<std::uint8_t> VTI::decode_b64(std::string_view b64) {
  std::string b64_cleaned;
  b64_cleaned.reserve(b64.size());
  for (char ch : b64) {
    if (std::isspace(static_cast<unsigned char>(ch)))
      continue;
    b64_cleaned.push_back(ch);
  }
  std::string_view b64_view = b64_cleaned;

  std::vector<std::uint8_t> out(boost::beast::detail::base64::decoded_size(b64_view.size()));
  std::size_t pos = 0;
  while (!b64_view.empty()) {
    while (!b64_view.empty() && b64_view.front() == '=') {
      b64_view.remove_prefix(1);
    }
    if (b64_view.empty()) break;

    std::size_t chunk_size = b64_view.size();
    std::size_t eq_pos = b64_view.find('=');
    // std::string_view::find returns npos if "=" is not found
    if (eq_pos != std::string_view::npos) {
      chunk_size = eq_pos;
      while (chunk_size < b64_view.size() && b64_view[chunk_size] == '=') {
        ++chunk_size;
      }
    }
    std::string_view chunk = b64_view.substr(0, chunk_size);
    std::span<std::uint8_t> out_chunk = std::span<std::uint8_t>(out).subspan(pos);
    auto [bytes_written, chars_consumed] = boost::beast::detail::base64::decode(out_chunk.data(), chunk.data(), chunk.size());
    if (chars_consumed == 0)
      IO::error(kVTKError, "base64 decode failed (invalid character or malformed input)");
    pos += bytes_written;
    b64_view.remove_prefix(chars_consumed);
  }
  out.resize(pos);
  return out;
}

template<class T>
std::vector<T> VTI::view(const std::vector<uint8_t>& raw) const {
  if (raw.size() % sizeof(T) != 0)
    IO::error(kVTKError, "size mismatch");
  std::vector<T> out(raw.size() / sizeof(T));
  std::memcpy(out.data(), raw.data(), raw.size());
  return out;
}

template<typename T>
void VTI::allocate_and_convert(const DecodedBuffer& d, CFI_cdesc_t* desc) {
  std::function<T*(std::size_t)> allocate_fortran_string = [&](std::size_t n) -> T* {
    if (desc->base_addr)
      CFI_deallocate(desc);
    CFI_index_t lb = 1, ext = static_cast<CFI_index_t>(n);
    if (CFI_allocate(desc, &lb, &ext, 0) != CFI_SUCCESS)
      throw std::runtime_error("CFI allocation failed");
    return static_cast<T*>(desc->base_addr);
  };

  const std::string& type = d.vtk_type;
  if (type == "Int32") {
    const std::vector<int32_t> src = view<int32_t>(d.raw_bytes);
    T* dst = allocate_fortran_string(src.size());
    std::transform(src.begin(), src.end(), dst, [](int32_t v){ return static_cast<T>(v); });
  } else if (type == "Int64") {
    const std::vector<int64_t> src = view<int64_t>(d.raw_bytes);
    T* dst = allocate_fortran_string(src.size());
    std::transform(src.begin(), src.end(), dst, [](int64_t v){ return static_cast<T>(v); });
  } else if (type == "Float32") {
    const std::vector<float> src = view<float>(d.raw_bytes);
    T* dst = allocate_fortran_string(src.size());
    std::transform(src.begin(), src.end(), dst, [](float v){ return static_cast<T>(v); });
  } else if (type == "Float64") {
    const std::vector<double> src = view<double>(d.raw_bytes);
    T* dst = allocate_fortran_string(src.size());
    if constexpr (std::is_same_v<T,double>)
      std::memcpy(dst, src.data(), src.size() * sizeof(double));   // 1-to-1 copy
    else
      std::transform(src.begin(), src.end(), dst, [](double v){ return static_cast<T>(v); });
  } else {
    IO::error(kVTKError, "unknown VTK type '" + type + '\'');
  }
}

extern "C" {
  gsl::owner<VTI*> VTI__new(const char* vti_path) {
    return new VTI(vti_path);
  }

  void VTI__delete(gsl::owner<VTI*> vti) {
    delete vti;
  }

  void C_VTI_readDatasetInt(VTI* vti, const char* name, CFI_cdesc_t* desc) {
    vti->read_dataset_int(name, desc);
  }

  void C_VTI_readDatasetReal(VTI* vti, const char* name, CFI_cdesc_t* desc) {
    vti->read_dataset_real(name, desc);
  }

  // https://clang.llvm.org/extra/clang-tidy/checks/bugprone/easily-swappable-parameters.html
  void C_VTI_readGeometry(VTI* vti,
                          int* cells_ptr,
                          // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                          double* geom_size_ptr,
                          double* origin_ptr,
                          CFI_cdesc_t* labels_desc) {
    vti->read_geometry(cells_ptr, geom_size_ptr, origin_ptr, labels_desc);
  }
}

#endif
