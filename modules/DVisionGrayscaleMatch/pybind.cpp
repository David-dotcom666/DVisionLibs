#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11/ndarray_converter.h"
#include "src/GrayscaleMatch.h"

namespace py = pybind11;

PYBIND11_MODULE(DVisionGrayscaleMatch, m) {
    NDArrayConverter::init_numpy();
    
    py::class_<GrayscaleMatch>(m, "GrayscaleMatch")
        .def(py::init<>())
        .def("setTempimage", &GrayscaleMatch::setTempimage, py::arg("_temp") = cv::Mat(), py::arg("useMask") = 0, py::arg("_mask") = cv::Mat())
        .def("setScore", &GrayscaleMatch::setScore)
        .def("setMaxtargs", &GrayscaleMatch::setMaxtargs)
        .def("setPyramidLayer", &GrayscaleMatch::setPyramidLayer)
        .def("setAngle", &GrayscaleMatch::setAngle)
        .def("setMaxOverlap", &GrayscaleMatch::setMaxOverlap)
        .def("setPolarity", &GrayscaleMatch::setPolarity)
        .def("setSubpixel", &GrayscaleMatch::setSubpixel)
        .def("setMean", &GrayscaleMatch::setMean)
        .def("setSDV", &GrayscaleMatch::setSDV)
        .def("match", &GrayscaleMatch::match)
        .def("getResult", &GrayscaleMatch::getResult)
        ;

    m.doc() = "DVision Grayscale match, version: 1.0.0, time: 20240210";
}

