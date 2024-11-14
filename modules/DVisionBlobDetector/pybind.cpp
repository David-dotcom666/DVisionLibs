#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11/ndarray_converter.h"
#include "src/BlobDetector.h"

namespace py = pybind11;

PYBIND11_MODULE(DVisionBlobDetector, m) {
    NDArrayConverter::init_numpy();
    
    py::class_<BlobDetector>(m,"BlobDetector")
        .def(py::init<>())
        .def("setColor", &BlobDetector::setColor)
        .def("setmaxTargets", &BlobDetector::setmaxTargets)
        .def("setSort", &BlobDetector::setSort, py::arg("sortMethod") = 1, py::arg("sortFeature") = 0)
        .def("setThreshold", &BlobDetector::setThreshold, py::arg("autoThreshold") = 1, py::arg("minThreshold") = 155, py::arg("maxThreshold") = 255)
        .def("setArea", &BlobDetector::setArea, py::arg("filterByArea") = 1, py::arg("minArea") = 25, py::arg("maxArea") = INT_MAX)
        .def("setCircularity", &BlobDetector::setCircularity, py::arg("filterByCircularity") = 0, py::arg("minCircularity") = 0.1, py::arg("maxCircularity") = 1)
        .def("setInertiaRatio", &BlobDetector::setInertiaRatio, py::arg("filterByInertia") = 0, py::arg("minInertiaRatio") = 0.1, py::arg("maxInertiaRatio") = 1)
        .def("setRectangularity", &BlobDetector::setRectangularity, py::arg("filterByRectangularity") = 0, py::arg("minRectangularity") = 0.1, py::arg("maxRectangularity") = 1)
        .def("setPerimeter", &BlobDetector::setPerimeter, py::arg("filterByPerimeter") = 1, py::arg("minPerimeter") = 25, py::arg("maxPerimeter") = INT_MAX)
        .def("setCentroidOffset", &BlobDetector::setCentroidOffset, py::arg("filterByCentroidOffset") = 1, py::arg("minCentroidOffset") = 25, py::arg("maxCentroidOffset") = INT_MAX)
        .def("detect", &BlobDetector::detect, py::arg("_image") = cv::Mat(), py::arg("_mask") = cv::Mat())
        .def("getResults", &BlobDetector::getResults)
        .def("getContours", &BlobDetector::getContours)
        .def("getExternalContours", &BlobDetector::getExternalContours)
        .def("getMinRects", &BlobDetector::getMinRects)
        .def("getBinaryimage", &BlobDetector::getBinaryimage)
        .def("getResultBinary", &BlobDetector::getResultBinary)
        ;

    m.doc() = "DVision Blob Detector, version: 1.0.0, time: 20241114";
}

