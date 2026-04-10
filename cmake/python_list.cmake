set(FILE_HEADER
    "#include <pybind11/pybind11.h>
    #include <vector>
    #include <functional>
    namespace py = pybind11;")

set(FUNC_DECL "")
set(VARIABLE_DESC "{")

foreach(file ${PGO_PYTHON_SOURCES})
    set(FUNC_DECL "${FUNC_DECL}
      void ${file}_init(py::module &m);")

    set(VARIABLE_DESC "${VARIABLE_DESC}
      ${file}_init,")
endforeach()

set(VARIABLE_DESC "${VARIABLE_DESC} };")

file(WRITE ${PGO_PYLIST_FILE}
    "${FILE_HEADER}
    ${FUNC_DECL}
    std::vector<std::function<void(py::module &)>> initFunctions = ${VARIABLE_DESC}")