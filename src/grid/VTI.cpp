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
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <span>
#include <sstream>
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
#include "../utils.h"
#include "../IO.h"

constexpr int VTK_ERROR = 844;

namespace {

/**
 * @brief Fetch an XML attribute, return an empty string if it doesn't exist.
 *
 * @param node XML node.
 * @param key  Attribute name.
 * @return     A string with the attribute if it exists, otherwise an empty one.
 */
std::string get_attr(const pt::ptree& node, const char* key) {
  if (auto a = node.get_child_optional("<xmlattr>")) {
    if (auto v = a->get_optional<std::string>(key))
      return *v;
  }
  return {};
}

/**
 * @brief Validate VTI file format.
 *
 * @param root The VTKFile root element of the property tree.
 * @throws std::runtime_error if validation fails.
 */
void check_file_format(const pt::ptree& root) {
  const std::string type = get_attr(root, "type");
  if (type != "ImageData")
    IO::error(VTK_ERROR, "type is not ImageData (got '" + type + "')");

  const std::string byte_order = get_attr(root, "byte_order");
  if (byte_order.empty() || byte_order != "LittleEndian")
    IO::error(VTK_ERROR, "byte_order must be 'LittleEndian' (got '" + byte_order + "')");

  const std::string compressor = get_attr(root, "compressor");
  if (!compressor.empty() && compressor != "vtkZLibDataCompressor")
    IO::error(VTK_ERROR, "compressor is not vtkZLibDataCompressor (got '" + compressor + "')");
}

/**
 * @brief Convert the VTK byte stream into a vector of the target type.
 *
 * @tparam T  Target element type (int32_t, int64_t, or double).
 * @param d   Decoded VTK Dataarray plus type tag.
 * @return    Converted values.
 */
template <typename T>
std::vector<T> cast_buffer(const DecodedBuffer& d) {
  auto convert = [&]<typename SrcT>() -> std::vector<T> {
    if constexpr (std::is_floating_point_v<SrcT> && std::is_integral_v<T>)
      IO::error(VTK_ERROR, "cannot cast floating-point to integer");
    if (d.raw_bytes.size() % sizeof(SrcT) != 0)
      IO::error(VTK_ERROR, "size mismatch");
    const std::size_t n = d.raw_bytes.size() / sizeof(SrcT);
    std::vector<T> dst(n);
    if constexpr (std::is_same_v<SrcT, T>) {
      std::memcpy(dst.data(), d.raw_bytes.data(), d.raw_bytes.size());
    } else {
      std::vector<SrcT> src(n);
      std::memcpy(src.data(), d.raw_bytes.data(), d.raw_bytes.size());
      std::transform(src.begin(), src.end(), dst.begin(), [](SrcT v) {
        return static_cast<T>(v);
      });
    }
    return dst;
  };

  if (d.vtk_type == "Int32")
    return convert.template operator()<int32_t>();
  if (d.vtk_type == "Int64")
    return convert.template operator()<int64_t>();
  if (d.vtk_type == "Float32")
    return convert.template operator()<float>();
  if (d.vtk_type == "Float64")
    return convert.template operator()<double>();

  IO::error(VTK_ERROR, "unknown VTK type '" + d.vtk_type + '\'');
  return {};
}

} // namespace

VTI::VTI(const fs::path& file_path) {
  if (file_path.empty())
    IO::error(VTK_ERROR, "no valid geometry file path supplied");
  std::ifstream f(file_path, std::ios::binary);
  if (!f)
    IO::error(VTK_ERROR, std::string("cannot open file '") + file_path.string() + '\'');
  pt::read_xml(f, tree);
  auto root = tree.get_child_optional("VTKFile");
  if (!root)
    IO::error(VTK_ERROR, "missing <VTKFile> element");
  check_file_format(*root);
  auto image_data = root->get_child_optional("ImageData");
  if (!image_data)
    IO::error(VTK_ERROR, "missing <ImageData> element");
}

std::uint64_t VTI::read_word(const std::span<const uint8_t> word) {
  if (word.size() != sizeof(std::uint32_t) && word.size() != sizeof(std::uint64_t))
    IO::error(VTK_ERROR, "unsupported VTI word size");
  std::uint64_t w = 0;
  std::memcpy(&w, word.data(), word.size());
  return w;
}

std::vector<std::uint8_t> VTI::decode_compressed(const std::string_view b64_string, const std::size_t n_bytes_per_word) {
  const std::vector<std::uint8_t> decoded_vec = VTI::decode_b64(b64_string);
  const std::span<const std::uint8_t> decoded(decoded_vec);

  const std::size_t n_header_bytes = 3 * n_bytes_per_word; // VTK generates headers of size 3
  if (decoded.size() < n_header_bytes)
    IO::error(VTK_ERROR, "header for compressed VTI too short");
  std::span<const std::uint8_t> header = decoded.first(n_header_bytes);
  std::uint64_t n_blocks = VTI::read_word(header.subspan(0 * n_bytes_per_word, n_bytes_per_word));
  std::uint64_t block_uncompressed = VTI::read_word(header.subspan(1 * n_bytes_per_word, n_bytes_per_word));
  std::uint64_t last_block_size = VTI::read_word(header.subspan(2 * n_bytes_per_word, n_bytes_per_word));

  const std::size_t header_size = n_header_bytes + static_cast<std::size_t>(n_blocks) * n_bytes_per_word;
  if (decoded.size() < header_size)
    IO::error(VTK_ERROR, "missing header for compressed VTI");
  std::span<const std::uint8_t> c_table = decoded.subspan(n_header_bytes, header_size - n_header_bytes);
  std::vector<std::size_t> block_sizes(n_blocks);
  for (std::size_t i = 0; i < n_blocks; ++i) {
    std::uint64_t word = VTI::read_word(c_table.subspan(i * n_bytes_per_word, n_bytes_per_word));
    block_sizes.at(i) = static_cast<std::size_t>(word);
  }

  std::span<const std::uint8_t> deflated = decoded.subspan(header_size);
  const std::size_t total_uncompressed =
      (n_blocks > 1 ? (n_blocks - 1) * static_cast<std::size_t>(block_uncompressed) : 0) +
      (last_block_size ? static_cast<std::size_t>(last_block_size) : static_cast<std::size_t>(block_uncompressed));
  std::vector<std::uint8_t> res_vec(total_uncompressed);
  std::span<std::uint8_t> res(res_vec);
  std::size_t src_offset = 0;
  std::size_t dst_offset = 0;
  for (std::size_t block_idx = 0; block_idx < n_blocks; ++block_idx) {
    if (src_offset + block_sizes.at(block_idx) > deflated.size())
      IO::error(VTK_ERROR, "invalid VTI, defined size overflowing for block " + std::to_string(block_idx));

    const std::size_t uncompressed_len = (block_idx + 1 == n_blocks && last_block_size > 0)
                                             ? static_cast<std::size_t>(last_block_size)
                                             : static_cast<std::size_t>(block_uncompressed);
    uLongf dst = static_cast<uLongf>(uncompressed_len);
    int z = uncompress(res.subspan(dst_offset).data(),
                       &dst,
                       deflated.subspan(src_offset).data(),
                       static_cast<uLongf>(block_sizes.at(block_idx)));
    if (z != Z_OK || dst != static_cast<uLongf>(uncompressed_len))
      IO::error(VTK_ERROR, "zlib inflate failed on block " + std::to_string(block_idx));
    src_offset += block_sizes.at(block_idx);
    dst_offset += uncompressed_len;
  }
  return res_vec;
}

std::vector<std::uint8_t> VTI::decode_uncompressed(const std::string_view b64_string, const std::size_t n_bytes_per_word) {
  const std::vector<std::uint8_t> decoded_vec = VTI::decode_b64(b64_string);
  const std::span<const std::uint8_t> decoded(decoded_vec);

  std::vector<std::uint8_t> res_vec;
  std::size_t offset = 0;

  while (offset + n_bytes_per_word <= decoded.size()) {
    const std::uint64_t n_bytes = VTI::read_word(decoded.subspan(offset, n_bytes_per_word));
    offset += n_bytes_per_word;

    if (n_bytes == 0)
      break;
    if (offset + n_bytes > decoded.size())
      IO::error(VTK_ERROR, "VTI uncompressed: data block exceeds payload");
    std::span<const std::uint8_t> payload = decoded.subspan(offset, static_cast<std::size_t>(n_bytes));
    res_vec.insert(res_vec.end(), payload.begin(), payload.end());
    offset += static_cast<std::size_t>(n_bytes);
  }
  return res_vec;
}

template <typename T>
std::vector<T> VTI::read_dataset(const std::string_view label) const {
  return cast_buffer<T>(read_dataset_raw(label));
}

template std::vector<std::int32_t> VTI::read_dataset<std::int32_t>(std::string_view) const;
template std::vector<std::int64_t> VTI::read_dataset<std::int64_t>(std::string_view) const;
template std::vector<double> VTI::read_dataset<double>(std::string_view) const;

void VTI::read_geometry(std::span<int, VTI::DIM> cells,
                        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        std::span<double, VTI::DIM> geom_size,
                        std::span<double, VTI::DIM> origin,
                        std::vector<std::string>& labels) const {

  auto parse_ints = [](const std::string& s) {
    std::istringstream is(s);
    std::vector<int> out;
    int v = 0;
    while (is >> v)
      out.push_back(v);
    return out;
  };

  auto parse_3_doubles = [](const std::string& field_name, const std::string& s) {
    std::istringstream is(s);
    std::array<double, 3> out{};
    is >> out[0] >> out[1] >> out[2];
    if (!is)
      IO::error(VTK_ERROR, "bad numeric field for '" + field_name + "' (got '" + s + "')");
    return out;
  };

  auto root = tree.get_child_optional("VTKFile");
  auto img = root->get_child_optional("ImageData");
  std::string dir = get_attr(*img, "Direction");
  if (!dir.empty() && dir != "1 0 0 0 1 0 0 0 1")
    IO::error(VTK_ERROR, "unsupported 'Direction' (got '" + dir + "')");
  const std::string extent_str = get_attr(*img, "WholeExtent");
  if (extent_str.empty())
    IO::error(VTK_ERROR, "missing 'WholeExtent'");
  const std::vector<int> extent = parse_ints(extent_str);
  if (extent.size() != 2 * VTI::DIM || (extent.at(0) != 0 || extent.at(2) != 0 || extent.at(4) != 0))
    IO::error(VTK_ERROR, "invalid 'WholeExtent' (got '" + extent_str + "')");

  const std::array<double, VTI::DIM> spacing = parse_3_doubles("ImageData@Spacing", get_attr(*img, "Spacing"));
  const std::array<double, VTI::DIM> origin_vec = parse_3_doubles("ImageData@Origin", get_attr(*img, "Origin"));

  /* modern form for component-wise assignment below (does not work on macOS)
  auto cells_vec = extent | std::views::drop(1) | std::views::stride(2);
  std::copy(cells_vec.begin(), cells_vec.end(), cells.begin());
  */
  cells[0] = extent.at(1);
  cells[1] = extent.at(3);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
  cells[2] = extent.at(5);
  geom_size[0] = spacing[0] * cells[0];
  geom_size[1] = spacing[1] * cells[1];
  geom_size[2] = spacing[2] * cells[2];
  std::copy(origin_vec.begin(), origin_vec.end(), origin.begin());

  if (geom_size[0] <= 0 || geom_size[1] <= 0 || geom_size[2] <= 0)
    IO::error(VTK_ERROR, "one or more entries <= 0 for 'size'");
  if (cells[0] < 1 || cells[1] < 1 || cells[2] < 1)
    IO::error(VTK_ERROR, "one or more entries < 1 for 'cells'");

  if (auto cell_data = img->get_child_optional("Piece.CellData")) {
    labels.reserve(cell_data->size());
    for (const auto& child : *cell_data) {
      if (child.first != "DataArray")
        continue;
      std::string name = get_attr(child.second, "Name");
      if (name.empty())
        continue;
      if (std::find(labels.begin(), labels.end(), name) != labels.end()) {
        IO::error(VTK_ERROR, "repeated label '" + name + '\'');
      }
      labels.push_back(std::move(name));
    }
  }
}

DecodedBuffer VTI::read_dataset_raw(const std::string_view label) const {
  auto root = tree.get_child_optional("VTKFile");
  auto image_data = root->get_child_optional("ImageData");
  auto cell_data = image_data->get_child_optional("Piece.CellData");
  if (!cell_data)
    IO::error(VTK_ERROR, "missing <CellData> element");
  boost::optional<const pt::ptree&> data_array_node;
  for (const auto& child : *cell_data) {
    if (child.first == "DataArray") {
      auto name_attr = get_attr(child.second, "Name");
      if (name_attr == label) {
        data_array_node = child.second;
        break;
      }
    }
  }
  if (!data_array_node)
    IO::error(VTK_ERROR, "no DataArray with Name='" + std::string(label) + "' found");
  if (get_attr(*data_array_node, "format") != "binary")
    IO::error(VTK_ERROR, "DataArray '" + std::string(label) + "' is not binary");
  const std::string vtk_type = get_attr(*data_array_node, "type");
  if (vtk_type.empty())
    IO::error(VTK_ERROR, "DataArray missing 'type' attribute");
  const bool compressed = (get_attr(*root, "compressor") == "vtkZLibDataCompressor");
  const std::size_t n_bytes_per_word = (get_attr(*root, "header_type") == "UInt64") ? 8 : 4;

  std::string text_content = data_array_node->get_value<std::string>();
  std::vector<std::uint8_t> raw_bytes = compressed ? decode_compressed(text_content, n_bytes_per_word)
                                                   : decode_uncompressed(text_content, n_bytes_per_word);

  return DecodedBuffer{std::move(vtk_type), std::move(raw_bytes)};
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
    if (b64_view.empty())
      break;

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
      IO::error(VTK_ERROR, "base64 decode failed (invalid character or malformed input)");
    pos += bytes_written;
    b64_view.remove_prefix(chars_consumed);
  }
  out.resize(pos);
  return out;
}

extern "C" {
VTI* C_VTI_new(const CFI_cdesc_t* vti_path) {
  return new VTI(fs::path(utils::view_on_fortran_string(vti_path))); // NOLINT(cppcoreguidelines-owning-memory)
}

void C_VTI_delete(VTI* vti) {
  delete vti; // NOLINT(cppcoreguidelines-owning-memory)
}

void C_VTI_readDatasetInt(const VTI* vti, const CFI_cdesc_t* label, CFI_cdesc_t* array_out) {
  const auto sv = utils::view_on_fortran_string(label);
  if (array_out->type == CFI_type_int32_t) {
    utils::fortran_array_from_vector(vti->read_dataset<std::int32_t>(sv), array_out);
  } else if (array_out->type == CFI_type_int64_t) {
    utils::fortran_array_from_vector(vti->read_dataset<std::int64_t>(sv), array_out);
  } else {
    IO::error(VTK_ERROR, "unsupported integer type for dataset '" + std::string(sv) + "'");
  }
}

void C_VTI_readDatasetReal(const VTI* vti, const CFI_cdesc_t* label, CFI_cdesc_t* array_out) {
  utils::fortran_array_from_vector(vti->read_dataset<double>(utils::view_on_fortran_string(label)), array_out);
}

// https://clang.llvm.org/extra/clang-tidy/checks/bugprone/easily-swappable-parameters.html
void C_VTI_readGeometry(const VTI* vti,
                        int* cells_ptr,
                        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        double* geom_size_ptr,
                        double* origin_ptr,
                        CFI_cdesc_t* array_out) {
  std::vector<std::string> labels;
  vti->read_geometry(std::span<int, VTI::DIM>(cells_ptr, VTI::DIM),
                     std::span<double, VTI::DIM>(geom_size_ptr, VTI::DIM),
                     std::span<double, VTI::DIM>(origin_ptr, VTI::DIM),
                     labels);
  if (array_out != nullptr) {
    utils::fortran_array_from_vector_str(labels, array_out);
  }
}
}

#endif
