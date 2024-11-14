#include<opencv2/opencv.hpp>
#include<ReadBarcode.h>
#include"src/BarcodeRecognize.h"

using namespace cv;
using namespace std;

//inline ZXing::ImageView ImageViewFromMat(const cv::Mat& image)
//{
//	using ZXing::ImageFormat;
//
//	auto fmt = ImageFormat::None;
//	switch (image.channels()) {
//	case 1: fmt = ImageFormat::Lum; break;
//	case 3: fmt = ImageFormat::BGR; break;
//	case 4: fmt = ImageFormat::BGRX; break;
//	}
//
//	if (image.depth() != CV_8U || fmt == ImageFormat::None)
//		return { nullptr, 0, 0, ImageFormat::None };
//
//	return { image.data, image.cols, image.rows, fmt };
//}
//
//inline ZXing::Results ReadBarcodes(const cv::Mat& image, const ZXing::DecodeHints& hints = {})
//{
//	return ZXing::ReadBarcodes(ImageViewFromMat(image), hints);
//}

inline void DrawResult(cv::Mat& img, ZXing::Result res)
{
	auto pos = res.position();
	auto zx2cv = [](ZXing::PointI p) { return cv::Point(p.x, p.y); };
	auto contour = std::vector<cv::Point>{ zx2cv(pos[0]), zx2cv(pos[1]), zx2cv(pos[2]), zx2cv(pos[3]) };
	const auto* pts = contour.data();
	int npts = contour.size();

	cv::polylines(img, &pts, &npts, 1, true, CV_RGB(0, 255, 0));
	cv::putText(img, res.text(), zx2cv(pos[3]) + cv::Point(0, 20), cv::FONT_HERSHEY_DUPLEX, 0.5, CV_RGB(0, 255, 0));
}

int main(int argc, char* argv[])
{
	Mat src1 = imread("img/DVisionBarcode/bar2.bmp");
	imshow("Display window1", src1);
	BarcodeRecognize t1 = BarcodeRecognize();
	if (1)
	{
		t1.recognizes(src1);
		auto results = t1.getResults();
		cout << results.size() << endl;
		DrawResult(src1, results[1]);

	}
	else
	{	
		t1.recognize(src1);
		auto result = t1.getResult();
		cout << result.value().text() << endl;
		DrawResult(src1, result.value());
	}

	imshow("Display window", src1);
	
	waitKey();
	return 0;
}



