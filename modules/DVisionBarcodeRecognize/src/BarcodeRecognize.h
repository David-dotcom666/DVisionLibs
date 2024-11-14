#pragma once

#include "BarcodeFormat.h"

// Reader
#include "ReadBarcode.h"
#include "ZXAlgorithms.h"

// Writer
#include "BitMatrix.h"
#include "MultiFormatWriter.h"

#include <vector>
#include<opencv2/opencv.hpp>
#include <optional>

struct Params
{
	bool is_pure = false; // 如果输入仅包含完美对齐的条形码（生成的图像），则设置为 True,在这种情况下加快检测速度。
	bool try_rotate = true; // True:解码器在任何方向搜索条形码; False:不会搜索 90° / 270° 旋转的条形码。
	bool try_downscale = true; // True:解码器会扫描输入的缩小版本；False:只会在提供的分辨率中搜索
	ZXing::BarcodeFormats formats = ZXing::BarcodeFormats{}; // 要解码的格式。 如果“无”，则解码所有格式。
	ZXing::TextMode text_mode = ZXing::TextMode::HRI; // 指定控制原始字节内容如何转码为结果中的文本的 TextMode
	ZXing::Binarizer binarizer = ZXing::Binarizer::LocalAverage; // 二值化器用于在解码条形码之前转换图像
	ZXing::EanAddOnSymbol ean_add_on_symbol = ZXing::EanAddOnSymbol::Ignore; // 指定在扫描 EAN/UPC 代码时是忽略、读取还是需要 EAN-2/5 附加符号。
};

class BarcodeRecognize
{
public:
	BarcodeRecognize() {};
	~BarcodeRecognize() {};

public:
	void setFormats(ZXing::BarcodeFormats formats);
	void setTextMode(ZXing::TextMode text_mode);
	void setBinarizer(ZXing::Binarizer binarizer);
	void setEanAddOnSymbol(ZXing::EanAddOnSymbol ean_add_on_symbol);
	void setRotate(bool try_rotate);
	void setDownscale(bool try_downscale);
	void setPure(bool is_pure);

public:
	void recognize(cv::Mat _image);
	std::optional<ZXing::Result> getResult() { return m_result; }

	void recognizes(cv::Mat _image);
	ZXing::Results getResults() { return m_results; }


	cv::Mat write_barcode(ZXing::BarcodeFormat format, std::string text, int width, int height, int quiet_zone, int ec_level);

private:
	Params m_params = Params();
	std::optional<ZXing::Result> m_result;
	ZXing::Results m_results;
};

