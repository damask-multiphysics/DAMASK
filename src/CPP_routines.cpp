#include "ISO_Fortran_binding.h"
#ifdef BOOST
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/lexical_cast.hpp>

extern "C"{
  void get_uuid_cpp(CFI_cdesc_t *uuid, int *stat){
    const std::string uuid_tmp = boost::lexical_cast<std::string>(boost::uuids::random_generator()());
    const char *uuid_c = uuid_tmp.c_str();
    if (CFI_allocate(uuid, (CFI_index_t *)0, (CFI_index_t *)0, strlen(uuid_c)) == 0){
      memcpy(uuid->base_addr, uuid_c, strlen(uuid_c));
      *stat = 0;
    }
    else *stat = 1;
  }
}
#endif

