#pragma once

#include<opencv2/opencv.hpp>

struct Params
{
	Params() {};
	// 二值化参数
	int autoThreshold = 1;
	int minThreshold = 150;
	int maxThreshold = 255;

	// 连通域颜色 黑0/白1
	bool blobColor = 1;

	// 面积筛选
	bool filterByArea = true;
	int minArea = 25;
	int maxArea = INT_MAX;

	// 长短轴比筛选
	bool filterByInertia = true;
	float minInertiaRatio = 0.1f;
	float maxInertiaRatio = 1;

	// 凸性筛选 凸形状指的是该形状的任意两点之间的线段都完全位于该形状内部
	/*bool filterByConvexity = true;
	float minConvexity = 0.95f;
	float maxConvexity = std::numeric_limits<float>::max();*/

	// 矩形度筛选（面积与最小外接矩形之比）
	bool filterByRectangularity = false;
	float minRectangularity = 0.1f;
	float maxRectangularity = 1;

	// 周长筛选
	bool filterByPerimeter = false;
	int minPerimeter = 25;
	int maxPerimeter = INT_MAX;

	// 圆形度筛选
	bool filterByCircularity = false;
	float minCircularity = 0.1;
	float maxCircularity = 1;

	// 质心偏移距筛选
	bool filterByCentroidOffset = false;
	int minCentroidOffset = 25;
	int maxCentroidOffset = INT_MAX;

	// 查找个数，排序特征，排序方式
	int maxTargets = 0;
	int sortMethod = 1;  // 0：不排序，1：降序，2：升序
	int sortFeature = 0; // 0：面积，1：周长，2：圆形度，3：矩形度，4：轴比，5：外接矩形中心x，6：中心y

};


struct Result
{
	int label;										// 连通域标签
	cv::Point2f centroid;							// 质心
	int area;										// 面积
	int perimeter;					    			// 周长
	float circularity;								// 圆形度
	float rectangularity;							// 矩形度
	float axleratio;								// 长短轴之比
	int centroidOffset;			    				// 质心偏移距离
	cv::RotatedRect minRect;						// 最小外接矩形
	std::vector<cv::Point> externalContour;			// 外轮廓
	std::vector<std::vector<cv::Point>> Contour;    // 轮廓
};


class BlobDetector
{
public:
	BlobDetector() {}
    ~BlobDetector() {}

public:
	void setColor(bool blobColor);
	void setmaxTargets(int maxTargets);
	void setSort(int sortMethod, int sortFeature = 0);
	void setThreshold(int autoThreshold, int minThreshold = 150, int maxThreshold = 255);
	void setArea(bool filterByArea, int minArea = 25, int maxArea = INT_MAX);
	void setCircularity(bool filterByCircularity, float minCircularity = 0.1f, float maxCircularity = 1);
	void setInertiaRatio(bool filterByInertia, float minInertiaRatio = 0.1f, float maxInertiaRatio = 1);
	void setRectangularity(bool filterByRectangularity, float minRectangularity = 0.1f, float maxRectangularity = 1);
	void setPerimeter(bool filterByPerimeter, int minPerimeter = 25, int maxPerimeter = INT_MAX);
	void setCentroidOffset(bool filterByCentroidOffset, int minCentroidOffset = 25, int maxCentroidOffset = INT_MAX);

public:
	void detect(cv::Mat _image, cv::Mat _mask = cv::Mat());
	cv::Mat getResults();
	cv::Mat getContours();
	cv::Mat getExternalContours();
	cv::Mat getMinRects();
	cv::Mat getBinaryimage() {return m_binaryImage;}
	cv::Mat getResultBinary() { return m_resultBinary; }

private:
	Params m_params = Params();
	std::vector<Result> m_results;
	cv::Mat m_binaryImage;
	cv::Mat m_resultBinary;
	bool m_show = 0;
};

