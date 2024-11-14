#include "BarcodeRecognize.h"

using namespace std;
using namespace cv;
using namespace ZXing;

/*********************************************
*=================read=======================*
*********************************************/

void BarcodeRecognize::setFormats(BarcodeFormats formats)
{
	m_params.formats = formats;
}
void BarcodeRecognize::setTextMode(TextMode text_mode)
{
	m_params.text_mode = text_mode;
}
void BarcodeRecognize::setBinarizer(Binarizer binarizer)
{
	m_params.binarizer = binarizer;
}
void BarcodeRecognize::setEanAddOnSymbol(EanAddOnSymbol ean_add_on_symbol)
{
	m_params.ean_add_on_symbol = ean_add_on_symbol;
}
void BarcodeRecognize::setRotate(bool try_rotate)
{
	m_params.try_rotate = try_rotate;
}
void BarcodeRecognize::setDownscale(bool try_downscale)
{
	m_params.try_downscale = try_downscale;
}
void BarcodeRecognize::setPure(bool is_pure)
{
	m_params.is_pure = is_pure;
}

auto read_barcodes_impl(Mat _image, const BarcodeFormats& formats, bool try_rotate, bool try_downscale, TextMode text_mode,
	Binarizer binarizer, bool is_pure, EanAddOnSymbol ean_add_on_symbol, uint8_t max_number_of_symbols = 0xff)
{
	const auto hints = DecodeHints()
		.setFormats(formats)
		.setTryRotate(try_rotate)
		.setTryDownscale(try_downscale)
		.setTextMode(text_mode)
		.setBinarizer(binarizer)
		.setIsPure(is_pure)
		.setMaxNumberOfSymbols(max_number_of_symbols)
		.setEanAddOnSymbol(ean_add_on_symbol);

	Results a;
	Mat src = _image.clone();
	ImageFormat imgfmt = ImageFormat::None;
	if (src.channels() == 1)
	{
		imgfmt = ImageFormat::Lum;
	}
	else if (src.channels() == 3)
	{
		imgfmt = ImageFormat::RGB;
	}
	else if (src.channels() == 4)
	{
		imgfmt = ImageFormat::RGBX;
	}
	else
	{
		return a;
	}

	const auto height = src.rows;
	const auto width = src.cols;
	const auto channels = src.channels();
	if (imgfmt == ImageFormat::None) {
		if (channels == 1)
			imgfmt = ImageFormat::Lum;
		else if (channels == 3)
			imgfmt = ImageFormat::BGR;
		else
			return a;
	}

	const auto bytes = src.data;
	return ReadBarcodes({ bytes, width, height, imgfmt, width * channels, channels }, hints);
}


void BarcodeRecognize::recognize(cv::Mat _image)
{
	auto res = read_barcodes_impl(_image, m_params.formats, m_params.try_rotate, m_params.try_downscale, m_params.text_mode, m_params.binarizer, m_params.is_pure, m_params.ean_add_on_symbol, 1);
	m_result = res.empty() ? std::nullopt : std::optional(res.front());
}

void BarcodeRecognize::recognizes(cv::Mat _image)
{
	m_results = read_barcodes_impl(_image, m_params.formats, m_params.try_rotate, m_params.try_downscale, m_params.text_mode, m_params.binarizer, m_params.is_pure, m_params.ean_add_on_symbol);
}



/*********************************************
*=================write======================*
*********************************************/

cv::Mat BarcodeRecognize::write_barcode(BarcodeFormat format, std::string text, int width, int height, int quiet_zone, int ec_level)
{
	auto writer = MultiFormatWriter(format).setEncoding(CharacterSet::UTF8).setMargin(quiet_zone).setEccLevel(ec_level);
	auto bitmap = writer.encode(text, width, height);
	Mat result = Mat(bitmap.height(), bitmap.width(), CV_8UC1);
	for (size_t y = 0; y < result.rows; y++)
	{
		uchar* ptr = result.ptr(y);
		for (size_t x = 0; x < result.cols; x++)
		{
			ptr[x] = bitmap.get(narrow_cast<int>(x), narrow_cast<int>(y)) ? 0 : 255;
		}
	}
	return result;
}


