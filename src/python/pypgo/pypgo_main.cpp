#include "pgo_c.h"

#include <pybind11/pybind11.h>

#include <vector>
#include <functional>

class pypgoInit
{
public:
  pypgoInit()
  {
    pgo_init();
    // code initialization
  }
  ~pypgoInit()
  {
    // finalize
  }
};

namespace py = pybind11;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

extern std::vector<std::function<void(py::module &)>> initFunctions;

PYBIND11_MODULE(pypgo, m)
{
  static pypgoInit init;

  for (const auto &func : initFunctions) {
    func(m);
  }

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}