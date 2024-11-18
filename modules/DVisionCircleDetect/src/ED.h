#pragma once

#include <opencv2/opencv.hpp>

/// Special defines
#define EDGE_VERTICAL   1
#define EDGE_HORIZONTAL 2

#define ANCHOR_PIXEL  254
#define EDGE_PIXEL    255

#define LEFT  1
#define RIGHT 2
#define UP    3
#define DOWN  4

enum GradientOperator { PREWITT_OPERATOR = 101, SOBEL_OPERATOR = 102, SCHARR_OPERATOR = 103, LSD_OPERATOR = 104 };

struct StackNode {
	int r, c;   // starting pixel起始像素
	int parent; // parent chain (-1 if no parent)
	int dir;    // direction where you are supposed to go方向
};

// Used during Edge Linking在边缘链接期间使用
struct Chain {

	int dir;                   // Direction of the chain链的方向
	int len;                   // # of pixels in the chain
	int parent;                // Parent of this node (-1 if no parent)
	int children[2];           // Children of this node (-1 if no children)
	cv::Point* pixels;         // Pointer to the beginning of the pixels array指向像素数组开头的指针
};

class ED {

public://ED 函数重载，输入参数类型不同
	//输入图像，梯度算法，梯度阈值，锚点阈值，扫描间隔，最小路径长度，高斯滤波sigma，梯度幅值计算方法
	ED(cv::Mat _srcImage, GradientOperator _op = PREWITT_OPERATOR, int _gradThresh = 20, int _anchorThresh = 0, int _scanInterval = 1, int _minPathLen = 10, double _sigma = 1.0, bool _sumFlag = true);
	ED(short* gradImg, uchar* dirImg, int _width, int _height, int _gradThresh, int _anchorThresh, int _scanInterval = 1, int _minPathLen = 10, bool selectStableAnchors = true);
	ED();//上面都是构造函数

	std::vector<std::vector<cv::Point>> getSegments();//边缘段
	std::vector<std::vector<cv::Point>> getSortedSegments();//边缘段排序

protected:
	int width; // 图像宽度
	int height; // 图像高度
	uchar* srcImg;
	std::vector<std::vector< cv::Point> > segmentPoints;//边缘段点
	double sigma; // Gaussian sigma
	cv::Mat smoothImage;//平滑图像
	uchar* edgeImg; // pointer to edge image data
	uchar* smoothImg; // pointer to smoothed image data
	int segmentNos;//边缘段编号
	int minPathLen;//最小长度
	cv::Mat srcImage;

private:
	void ComputeGradient();//计算梯度
	void ComputeAnchorPoints();//提取锚点
	void JoinAnchorPointsUsingSortedAnchors();//连接锚点
	void sortAnchorsByGradValue();//按梯度值排序锚点
	int* sortAnchorsByGradValue1();//按梯度值排序锚点

	static int LongestChain(Chain* chains, int root);//最长链
	static int RetrieveChainNos(Chain* chains, int root, int chainNos[]);//边缘段总数

	int anchorNos;//锚点总数
	std::vector<cv::Point> anchorPoints;//锚点容器
	std::vector<cv::Point> edgePoints;//边缘点容器

	cv::Mat edgeImage;//边缘图像
	cv::Mat gradImage;//梯度图像

	uchar* dirImg; // 方向图像数据的指针
	short* gradImg; // 梯度图像数据的指针

	GradientOperator op; // 梯度计算方法（sobel等）
	int gradThresh; // 梯度阈值
	int anchorThresh; // 锚点阈值（影响锚点数量）
	int scanInterval;//扫描间隔（不超过锚点选择窗口的半径）
	bool sumFlag;
};

