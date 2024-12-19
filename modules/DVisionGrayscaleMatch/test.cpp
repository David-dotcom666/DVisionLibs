#include "src/GrayscaleMatch.h"

using namespace std;
using namespace cv;

int main()
{
	Mat src_image, temp_image, src_image2;
	if (1)
	{
		src_image = imread("img/DVisionGrayscaleMatch/1src_image/img1.bmp", 0);
		temp_image = imread("img/DVisionGrayscaleMatch/2temp_image/img1.bmp", 0);
		src_image2 = imread("img/DVisionGrayscaleMatch/1src_image/img1.bmp");
		Scalar meanMat, meanStdMat;
		meanStdDev(temp_image, meanMat, meanStdMat);
	}
	else
	{
		src_image = imread("img/DVisionGrayscaleMatch/image/6_1.bmp", 0);
		temp_image = imread("img/DVisionGrayscaleMatch/image/6_1.png", 0);
		src_image2 = imread("img/DVisionGrayscaleMatch/image/6_1.bmp");
	}

	double t1 = getTickCount();
	GrayscaleMatch m_GrayscaleMatch;
	m_GrayscaleMatch.setTempimage(temp_image);
	m_GrayscaleMatch.setScore(0.8);
	m_GrayscaleMatch.setAngle(0, 360);
	m_GrayscaleMatch.match(src_image2);
	cv::Mat result = m_GrayscaleMatch.getResult();
	double t2 = getTickCount();
	std::cout << "Match time: " << (t2 - t1) / getTickFrequency() << "s" << endl;

	std::cout << result.rows << endl;
	for (int i = 0; i < result.rows; i++)
	{
		int k = 2;
		float* rptr = (float*)result.ptr(i);
		circle(src_image2, Point((int)rptr[0], (int)rptr[1]), 1, { 0, 0, 255 }, 1);
		line(src_image2, Point((int)rptr[10], (int)rptr[11]), Point((int)rptr[4], (int)rptr[5]), Scalar(0, 255, 0), k, LINE_AA);
		line(src_image2, Point((int)rptr[4], (int)rptr[5]), Point((int)rptr[6], (int)rptr[7]), Scalar(0, 255, 0), k, LINE_AA);
		line(src_image2, Point((int)rptr[6], (int)rptr[7]), Point((int)rptr[8], (int)rptr[9]), Scalar(0, 255, 0), k, LINE_AA);
		line(src_image2, Point((int)rptr[8], (int)rptr[9]), Point((int)rptr[10], (int)rptr[11]), Scalar(0, 255, 0), k, LINE_AA);
	}
	imshow("src_image2", src_image2);
	//imwrite("src_image2.bmp", src_image2);

	cv::waitKey();
	return 0;
}


