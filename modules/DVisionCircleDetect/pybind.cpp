#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11/ndarray_converter.h"
#include"src/CircleDetect.h"


namespace py = pybind11;

PYBIND11_MODULE(DVisionCircleDetect, m) {
    NDArrayConverter::init_numpy();
    
    py::class_<CircleDetect>(m,"CircleDetect")
        .def(py::init<>())
        .def("setSplit", &CircleDetect::setSplit)
        .def("setInlier", &CircleDetect::setInlier)
        .def("setDistance", &CircleDetect::setDistance)
        .def("setFiltSize", &CircleDetect::setFiltSize)
        .def("setMinLen", &CircleDetect::setMinLen)
        .def("detect", &CircleDetect::detect)
        .def("getResult", &CircleDetect::getResult)
        ;

    m.doc() = "DVision Circle Detect, version: 1.0.0, time: 20230725";
}

