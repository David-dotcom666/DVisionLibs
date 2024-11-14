#include "src/BlobDetector.h"

using namespace cv;
using namespace std;

int main()
{

	int imgid = 3;
	Mat image2;
	switch (imgid)
	{
	case 1:
	{
		image2 = imread("img/DVisionBlobDetector/image/640.png", 0);
		break;
	}
	case 2:
	{
		image2 = imread("img/DVisionBlobDetector/image/a1.bmp", 0);
		break;
	}
	case 3:
	{
		image2 = imread("img/DVisionBlobDetector/image/2DMeasure.bmp");
		break;
	}
	}


	double t1 = getTickCount();
	BlobDetector m_BlobDetector = BlobDetector();
	// 阈值
	m_BlobDetector.setThreshold(1, 70, 255);
	// 连通域颜色
	m_BlobDetector.setColor(1);
	// 面积筛选
	m_BlobDetector.setArea(1, 2500);
	// 周长筛选
	m_BlobDetector.setPerimeter(1, 50);
	// 轴比筛选
	m_BlobDetector.setInertiaRatio(0, 0.1);
	// 矩形度筛选
	m_BlobDetector.setRectangularity(0, 0.9);
	// 圆形度筛选
	m_BlobDetector.setCircularity(0, 0.8);
	// 质心距离筛选
	m_BlobDetector.setCentroidOffset(0, 2);
	// 最大数量
	m_BlobDetector.setmaxTargets(0);
	// 排序
	m_BlobDetector.setSort(1,0);


	m_BlobDetector.detect(image2);
	// 获取连通域结果参数
	Mat result = m_BlobDetector.getResults();
	// 获取连通域最小外接矩形参数
	Mat minrect = m_BlobDetector.getMinRects();
	// 获取连通域所有轮廓
	Mat contour = m_BlobDetector.getContours();
	// 获取连通域最外层轮廓
	Mat externContour = m_BlobDetector.getExternalContours();

	cout << result.size() << endl;
	cout << minrect.size() << endl;
	cout << contour.size() << endl;
	cout << externContour.size() << endl;

	double t2 = getTickCount();
	std::cout << "BlobDetector time: " << (t2 - t1) / getTickFrequency() << "s" << endl;





	cv::waitKey();
	return 0;
}




