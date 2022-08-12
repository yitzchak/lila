#include <clasp/clasp.h>
#include <clasp/core/translators.h>
#include <lila/lila.h>

namespace lila {

/*CL_LAMBDA(core:&va-rest args)
CL_NAME("real-single-vector");
CL_DEFUN RealSingleVector_sp RealSingleVector_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<RealSingleVector_O>::allocate_with_default_constructor();
  res->_Value.resize(args->total_nargs());
  for (size_t i = 0; args->remaining_nargs() > 0; i++) {
    res->_Value[i] = core::clasp_to_float(args->next_arg());
  }
  return res;
}*/

} // Namespace lila
