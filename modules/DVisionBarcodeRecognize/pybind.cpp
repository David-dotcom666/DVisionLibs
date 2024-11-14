#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/ndarray_converter.h>
#include "src/BarcodeRecognize.h"

namespace py = pybind11;
using namespace cv;
using namespace ZXing;

std::ostream& operator<<(std::ostream& os, const Position& points) {
	for (const auto& p : points)
		os << p.x << "x" << p.y << " ";
	os.seekp(-1, os.cur);
	os << '\0';
	return os;
}


PYBIND11_MODULE(DVisionBarcodeRecognize, m) {
	NDArrayConverter::init_numpy();

	py::class_<BarcodeFormats> pyBarcodeFormats(m, "BarcodeFormats");

	py::enum_<BarcodeFormat>(m, "BarcodeFormat", py::arithmetic{}, "Enumeration of zxing supported barcode formats")
		.value("Aztec", BarcodeFormat::Aztec)
		.value("Codabar", BarcodeFormat::Codabar)
		.value("Code39", BarcodeFormat::Code39)
		.value("Code93", BarcodeFormat::Code93)
		.value("Code128", BarcodeFormat::Code128)
		.value("DataMatrix", BarcodeFormat::DataMatrix)
		.value("EAN8", BarcodeFormat::EAN8)
		.value("EAN13", BarcodeFormat::EAN13)
		.value("ITF", BarcodeFormat::ITF)
		.value("MaxiCode", BarcodeFormat::MaxiCode)
		.value("PDF417", BarcodeFormat::PDF417)
		.value("QRCode", BarcodeFormat::QRCode)
		.value("MicroQRCode", BarcodeFormat::MicroQRCode)
		.value("DataBar", BarcodeFormat::DataBar)
		.value("DataBarExpanded", BarcodeFormat::DataBarExpanded)
		.value("UPCA", BarcodeFormat::UPCA)
		.value("UPCE", BarcodeFormat::UPCE)
		.value("NONE", BarcodeFormat::None)
		.value("LinearCodes", BarcodeFormat::LinearCodes)
		.value("MatrixCodes", BarcodeFormat::MatrixCodes)
		.export_values()
		.def("__or__", [](BarcodeFormat f1, BarcodeFormat f2) { return f1 | f2; });
	pyBarcodeFormats
		.def("__repr__", py::overload_cast<BarcodeFormats>(static_cast<std::string(*)(BarcodeFormats)>(ToString)))
		.def("__str__", py::overload_cast<BarcodeFormats>(static_cast<std::string(*)(BarcodeFormats)>(ToString)))
		.def("__eq__", [](BarcodeFormats f1, BarcodeFormats f2) { return f1 == f2; })
		.def("__or__", [](BarcodeFormats fs, BarcodeFormat f) { return fs | f; })
		.def(py::init<BarcodeFormat>());
	py::implicitly_convertible<BarcodeFormat, BarcodeFormats>();

	py::enum_<Binarizer>(m, "Binarizer", "Enumeration of binarizers used before decoding images")
		.value("BoolCast", Binarizer::BoolCast)
		.value("FixedThreshold", Binarizer::FixedThreshold)
		.value("GlobalHistogram", Binarizer::GlobalHistogram)
		.value("LocalAverage", Binarizer::LocalAverage)
		.export_values();
	py::enum_<EanAddOnSymbol>(m, "EanAddOnSymbol", "Enumeration of options for EAN-2/5 add-on symbols check")
		.value("Ignore", EanAddOnSymbol::Ignore, "Ignore any Add-On symbol during read/scan")
		.value("Read", EanAddOnSymbol::Read, "Read EAN-2/EAN-5 Add-On symbol if found")
		.value("Require", EanAddOnSymbol::Require, "Require EAN-2/EAN-5 Add-On symbol to be present")
		.export_values();
	py::enum_<ContentType>(m, "ContentType", "Enumeration of content types")
		.value("Text", ContentType::Text)
		.value("Binary", ContentType::Binary)
		.value("Mixed", ContentType::Mixed)
		.value("GS1", ContentType::GS1)
		.value("ISO15434", ContentType::ISO15434)
		.value("UnknownECI", ContentType::UnknownECI)
		.export_values();
	py::enum_<TextMode>(m, "TextMode", "")
		.value("Plain", TextMode::Plain, "bytes() transcoded to unicode based on ECI info or guessed charset (the default mode prior to 2.0)")
		.value("ECI", TextMode::ECI, "standard content following the ECI protocol with every character set ECI segment transcoded to unicode")
		.value("HRI", TextMode::HRI, "Human Readable Interpretation (dependent on the ContentType)")
		.value("Hex", TextMode::Hex, "bytes() transcoded to ASCII string of HEX values")
		.value("Escaped", TextMode::Escaped, "Use the EscapeNonGraphical() function (e.g. ASCII 29 will be transcoded to '<GS>'")
		.export_values();
	py::class_<PointI>(m, "Point", "Represents the coordinates of a point in an image")
		.def_readonly("x", &PointI::x,
			":return: horizontal coordinate of the point\n"
			":rtype: int")
		.def_readonly("y", &PointI::y,
			":return: vertical coordinate of the point\n"
			":rtype: int");
	py::class_<Position>(m, "Position", "The position of a decoded symbol")
		.def_property_readonly("top_left", &Position::topLeft,
			":return: coordinate of the symbol's top-left corner\n"
			":rtype: zxing.Point")
		.def_property_readonly("top_right", &Position::topRight,
			":return: coordinate of the symbol's top-right corner\n"
			":rtype: zxing.Point")
		.def_property_readonly("bottom_left", &Position::bottomLeft,
			":return: coordinate of the symbol's bottom-left corner\n"
			":rtype: zxing.Point")
		.def_property_readonly("bottom_right", &Position::bottomRight,
			":return: coordinate of the symbol's bottom-right corner\n"
			":rtype: zxing.Point")
		.def("__str__", [](Position pos) {
		std::ostringstream oss;
		oss << pos;
		return oss.str();
			});
	py::class_<Result>(m, "Result", "Result of barcode reading")
		.def_property_readonly("valid", &Result::isValid,
			":return: whether or not result is valid (i.e. a symbol was found)\n"
			":rtype: bool")
		.def_property_readonly("text", [](const Result& res) { return res.text(); },
			":return: text of the decoded symbol (see also TextMode parameter)\n"
			":rtype: str")
		.def_property_readonly("bytes", [](const Result& res) { return py::bytes(res.bytes().asString()); },
			":return: uninterpreted bytes of the decoded symbol\n"
			":rtype: bytes")
		.def_property_readonly("format", &Result::format,
			":return: decoded symbol format\n"
			":rtype: zxing.BarcodeFormat")
		.def_property_readonly("symbology_identifier", &Result::symbologyIdentifier,
			":return: decoded symbology idendifier\n"
			":rtype: str")
		.def_property_readonly("content_type", &Result::contentType,
			":return: content type of symbol\n"
			":rtype: zxing.ContentType")
		.def_property_readonly("position", &Result::position,
			":return: position of the decoded symbol\n"
			":rtype: zxing.Position")
		.def_property_readonly("orientation", &Result::orientation,
			":return: orientation (in degree) of the decoded symbol\n"
			":rtype: int");
	m.def("barcode_format_from_str", &BarcodeFormatFromString,
		py::arg("str"),
		"Convert string to BarcodeFormat\n\n"
		":type str: str\n"
		":param str: string representing barcode format\n"
		":return: corresponding barcode format\n"
		":rtype: zxing.BarcodeFormat");
	m.def("barcode_formats_from_str", &BarcodeFormatsFromString,
		py::arg("str"),
		"Convert string to BarcodeFormats\n\n"
		":type str: str\n"
		":param str: string representing a list of barcodes formats\n"
		":return: corresponding barcode formats\n"
		":rtype: zxing.BarcodeFormats");

	py::class_<BarcodeRecognize>(m, "BarcodeRecognize")
		.def(py::init<>())
		.def("setFormats", &BarcodeRecognize::setFormats)
		.def("setTextMode", &BarcodeRecognize::setTextMode)
		.def("setBinarizer", &BarcodeRecognize::setBinarizer)
		.def("setEanAddOnSymbol", &BarcodeRecognize::setEanAddOnSymbol)
		.def("setRotate", &BarcodeRecognize::setRotate)
		.def("setDownscale", &BarcodeRecognize::setDownscale)
		.def("setPure", &BarcodeRecognize::setPure)
		.def("recognize", &BarcodeRecognize::recognize)
		.def("recognizes", &BarcodeRecognize::recognizes)
		.def("getResult", &BarcodeRecognize::getResult)
		.def("getResults", &BarcodeRecognize::getResults)
		;

	m.doc() = "DVision BarCode recognition, version: 1.0.0, time: 20241114";
}


