#include"src/CircleDetect.h"

using namespace cv;
using namespace std;

int main()
{
	Mat image = imread("img/DVisionShapeDetect/2.jpg");

	CircleDetect m_CircleDetect;
	m_CircleDetect.detect(image);
	cv::Mat result = m_CircleDetect.getResult();

	for (int i = 0; i < result.rows; i++)
	{
		float* IMG = (float*)result.ptr(i);
		circle(image, Point2d((int)IMG[0], (int)IMG[1]), (int)IMG[2], Scalar(0, 255, 0), 2);
	}
	imshow("result", image);
	imwrite("result.jpg", image);
	cv::waitKey();
	return 0;
}

