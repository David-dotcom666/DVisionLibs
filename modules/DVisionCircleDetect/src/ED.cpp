#include "ED.h"
#include <fstream>
#include"InstructionSet.h"

using namespace cv;
using namespace std;
//类ED的构造函数ED
//输入图像，梯度算法，梯度阈值，锚点阈值（梯度的差，影响锚点数量），扫描间隔，最小路径长度，高斯滤波sigma，梯度幅值计算方法
//一些参数，以及获取平滑图像，梯度图像，边缘方向图，边缘图
ED::ED(Mat _srcImage, GradientOperator _op, int _gradThresh, int _anchorThresh, int _scanInterval, int _minPathLen, double _sigma, bool _sumFlag){
	// Check parameters for sanity
	if (_gradThresh < 1) _gradThresh = 1;
	if (_anchorThresh < 0) _anchorThresh = 0;
	if (_sigma < 1.0) _sigma = 1.0;

	srcImage = _srcImage;
	height = srcImage.rows;
	width = srcImage.cols;
	op = _op;
	gradThresh = _gradThresh;
	anchorThresh = _anchorThresh;
	scanInterval = _scanInterval;
	minPathLen = _minPathLen;
	sigma = _sigma;
	sumFlag = _sumFlag;

	segmentNos = 0;
	segmentPoints.push_back(vector<Point>()); // create empty vector of points for segments

	edgeImage = Mat(height, width, CV_8UC1, Scalar(0)); // initialize edge Image//CV_8UC1:8bites Unsign C1:灰度图，Scalar(0):初始化值
	smoothImage = Mat(height, width, CV_8UC1);
	gradImage = Mat(height, width, CV_16SC1); // gradImage contains short values 

	srcImg = srcImage.data;

	//// Detect Edges By Edge Drawing Algorithm  ////

	/*------------ 高斯滤波平滑图像 -------------------*/
	if (sigma == 1.0)
		GaussianBlur(srcImage, smoothImage, Size(3, 3), sigma);
	else
		GaussianBlur(srcImage, smoothImage, Size(), sigma); // calculate kernel from sigma

	// Assign Pointers from Mat's data
	smoothImg = smoothImage.data;
	gradImg = (short*)gradImage.data;
	edgeImg = edgeImage.data;

	dirImg = new unsigned char[width * height];//边缘方向

	/*------------计算梯度幅值和边缘方向图-------------------*/
	ComputeGradient();

	/*------------ 提取锚点 -------------------*/
	ComputeAnchorPoints();

	/*------------连接锚点 -------------------*/
	JoinAnchorPointsUsingSortedAnchors();

	delete[] dirImg;
}

ED::ED(short* _gradImg, uchar* _dirImg, int _width, int _height, int _gradThresh, int _anchorThresh, int _scanInterval, int _minPathLen, bool selectStableAnchors)
{
	height = _height;
	width = _width;

	gradThresh = _gradThresh;
	anchorThresh = _anchorThresh;
	scanInterval = _scanInterval;
	minPathLen = _minPathLen;

	gradImg = _gradImg;
	dirImg = _dirImg;

	edgeImage = Mat(height, width, CV_8UC1, Scalar(0)); // initialize edge Image

	edgeImg = edgeImage.data;

	if (selectStableAnchors) {

		// Compute anchors with the user supplied parameters
		anchorThresh = 0; // anchorThresh used as zero while computing anchor points if selectStableAnchors set. 
						  // Finding higher number of anchors is OK, because we have following validation steps in selectStableAnchors.
		ComputeAnchorPoints();
		anchorThresh = _anchorThresh; // set it to its initial argument value for further anchor validation.
		anchorPoints.clear(); // considering validation step below, it should constructed again.

		for (int i = 1; i < height - 1; i++) {
			for (int j = 1; j < width - 1; j++) {
				if (edgeImg[i * width + j] != ANCHOR_PIXEL) continue;

				// Take only "stable" anchors
				// 0 degree edge
				if (edgeImg[i * width + j - 1] && edgeImg[i * width + j + 1]) {
					int diff1 = gradImg[i * width + j] - gradImg[(i - 1) * width + j];
					int diff2 = gradImg[i * width + j] - gradImg[(i + 1) * width + j];
					if (diff1 >= anchorThresh && diff2 >= anchorThresh) edgeImg[i * width + j] = 255;

					continue;
				} //end-if

				  // 90 degree edge
				if (edgeImg[(i - 1) * width + j] && edgeImg[(i + 1) * width + j]) {
					int diff1 = gradImg[i * width + j] - gradImg[i * width + j - 1];
					int diff2 = gradImg[i * width + j] - gradImg[i * width + j + 1];
					if (diff1 >= anchorThresh && diff2 >= anchorThresh) edgeImg[i * width + j] = 255;

					continue;
				} //end-if

				  // 135 degree diagonal
				if (edgeImg[(i - 1) * width + j - 1] && edgeImg[(i + 1) * width + j + 1]) {
					int diff1 = gradImg[i * width + j] - gradImg[(i - 1) * width + j + 1];
					int diff2 = gradImg[i * width + j] - gradImg[(i + 1) * width + j - 1];
					if (diff1 >= anchorThresh && diff2 >= anchorThresh) edgeImg[i * width + j] = 255;
					continue;
				} //end-if

				  // 45 degree diagonal
				if (edgeImg[(i - 1) * width + j + 1] && edgeImg[(i + 1) * width + j - 1]) {
					int diff1 = gradImg[i * width + j] - gradImg[(i - 1) * width + j - 1];
					int diff2 = gradImg[i * width + j] - gradImg[(i + 1) * width + j + 1];
					if (diff1 >= anchorThresh && diff2 >= anchorThresh) edgeImg[i * width + j] = 255;
				} //end-if

			} //end-for
		} //end-for

		for (int i = 0; i < width * height; i++)
			if (edgeImg[i] == ANCHOR_PIXEL)
				edgeImg[i] = 0;
			else if (edgeImg[i] == 255) {
				edgeImg[i] = ANCHOR_PIXEL;
				int y = i / width;
				int x = i % width;
				anchorPoints.push_back(Point(x, y)); // push validated anchor point to vector
			}

		anchorNos = anchorPoints.size(); // get # of anchor pixels
	}

	else {
		// Compute anchors with the user supplied parameters
		ComputeAnchorPoints(); // anchorThresh used as given as argument. No validation applied. (No stable anchors.)
	} //end-else

	segmentNos = 0;
	segmentPoints.push_back(vector<Point>()); // create empty vector of points for segments

	JoinAnchorPointsUsingSortedAnchors();
}

ED::ED()
{
	//
}

std::vector<std::vector<Point>> ED::getSegments()
{
	return segmentPoints;
}

//边缘段排序
std::vector<std::vector<Point>> ED::getSortedSegments()
{
	// sort segments from largest to smallest
	std::sort(segmentPoints.begin(), segmentPoints.end(), [](const std::vector<Point>& a, const std::vector<Point>& b) { return a.size() > b.size(); });

	return segmentPoints;
}

//计算梯度和边缘方向图，SSE版和原版

void ED::ComputeGradient()
{

	//Initialize gradient image for row = 0, row = height-1, column=0, column=width-1 
   //初始化梯度矩阵的第一行和最后一行
	for (int j = 0; j < width; j++) { gradImg[j] = gradImg[(height - 1) * width + j] = gradThresh - 1; }
	//初始化第一列和最后一列
	for (int i = 1; i < height - 1; i++) { gradImg[i * width] = gradImg[(i + 1) * width - 1] = gradThresh - 1; }
#pragma omp parallel for schedule(static)
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {
			int com1 = smoothImg[(i + 1) * width + j + 1] - smoothImg[(i - 1) * width + j - 1];
			int com2 = smoothImg[(i - 1) * width + j + 1] - smoothImg[(i + 1) * width + j - 1];
			int gx;
			int gy;
			//根据不同方法计算梯度图
			switch (op)
			{
			case PREWITT_OPERATOR:
				gx = abs(com1 + com2 + (smoothImg[i * width + j + 1] - smoothImg[i * width + j - 1]));
				gy = abs(com1 - com2 + (smoothImg[(i + 1) * width + j] - smoothImg[(i - 1) * width + j]));
				break;
			case SOBEL_OPERATOR:
				gx = abs(com1 + com2 + 2 * (smoothImg[i * width + j + 1] - smoothImg[i * width + j - 1]));
				gy = abs(com1 - com2 + 2 * (smoothImg[(i + 1) * width + j] - smoothImg[(i - 1) * width + j]));
				break;
			case SCHARR_OPERATOR:
				gx = abs(3 * (com1 + com2) + 10 * (smoothImg[i * width + j + 1] - smoothImg[i * width + j - 1]));
				gy = abs(3 * (com1 - com2) + 10 * (smoothImg[(i + 1) * width + j] - smoothImg[(i - 1) * width + j]));
				break;
			case LSD_OPERATOR:
				// com1 and com2 differs from previous operators, because LSD has 2x2 kernel
				int com1 = smoothImg[(i + 1) * width + j + 1] - smoothImg[i * width + j];
				int com2 = smoothImg[i * width + j + 1] - smoothImg[(i + 1) * width + j];
				gx = abs(com1 + com2);
				gy = abs(com1 - com2);
				break;
			}

			int sum;//梯度幅值
			//sumFlag:梯度幅值计算方法
			if (sumFlag)
				sum = gx + gy;
			else
				sum = (int)sqrt((double)gx * gx + (double)gy * gy);

			int index = i * width + j;//当前像素位置（data）
			gradImg[index] = sum;
			//梯度阈值与边缘方向
			if (sum >= gradThresh) {
				if (gx >= gy) dirImg[index] = EDGE_VERTICAL;//垂直方向
				else          dirImg[index] = EDGE_HORIZONTAL;//水平方向
			} //end-if
		} // end-for
	} // end-for
}

//计算梯度幅值与边缘方向
//void ComputeGradient()
//{
//#ifdef _SSE_
//	for (int i = 1; i < height - 1; i++) {
//		int j = 1;
//		for (; j < width - 16; j += 16) {
//			//前64位
//			__m128i p1 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j - 1))));
//			__m128i p2 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j))));
//			__m128i p3 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j + 1))));
//			__m128i p4 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + (i * width + j - 1))));
//			__m128i p6 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + (i * width + j + 1))));
//			__m128i p7 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j - 1))));
//			__m128i p8 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j))));
//			__m128i p9 = _mm_cvtepu8_epi16(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j + 1))));
//			//后64位
//			__m128i pb = { 8, 9, 10, 11, 12, 13, 14, 15 , 0, 1, 2, 3, 4, 5, 6, 7 };
//			__m128i p11 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j - 1))), pb));
//			__m128i p22 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j))), pb));
//			__m128i p33 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i - 1) * width + j + 1))), pb));
//			__m128i p44 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + (i * width + j - 1))), pb));
//			__m128i p66 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + (i * width + j + 1))), pb));
//			__m128i p77 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j - 1))), pb));
//			__m128i p88 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j))), pb));
//			__m128i p99 = _mm_cvtepu8_epi16(_mm_shuffle_epi8(_mm_loadu_si128((__m128i*) (smoothImg + ((i + 1) * width + j + 1))), pb));
//			//gx, gy
//			__m128i gx_front8 = _mm_abs_epi16(_mm_sub_epi16(_mm_add_epi16(_mm_add_epi16(p3, p6), p9), _mm_add_epi16(_mm_add_epi16(p1, p4), p7)));
//			__m128i gy_front8 = _mm_abs_epi16(_mm_sub_epi16(_mm_add_epi16(_mm_add_epi16(p7, p8), p9), _mm_add_epi16(_mm_add_epi16(p1, p2), p3)));
//			__m128i gx_back8 = _mm_abs_epi16(_mm_sub_epi16(_mm_add_epi16(_mm_add_epi16(p33, p66), p99), _mm_add_epi16(_mm_add_epi16(p11, p44), p77)));
//			__m128i gy_back8 = _mm_abs_epi16(_mm_sub_epi16(_mm_add_epi16(_mm_add_epi16(p77, p88), p99), _mm_add_epi16(_mm_add_epi16(p11, p22), p33)));
//			__m128i g_front8 = _mm_add_epi16(gx_front8, gy_front8);
//			__m128i g_back8 = _mm_add_epi16(gx_back8, gy_back8);
//			//方向
//			__m128i dre_front8 = _mm_cmplt_epi16(gx_front8, gy_front8);
//			__m128i dre_back8 = _mm_cmplt_epi16(gx_back8, gy_back8);
//			//结果
//			_mm_storeu_si128((__m128i*)(GradImg + (i * width + j)), g_front8);
//			_mm_storeu_si128((__m128i*)(GradImg + (i * width + j + 8)), g_back8);
//			_mm_storeu_si128((__m128i*)(DireImg + (i * width + j)), dre_front8);
//			_mm_storeu_si128((__m128i*)(DireImg + (i * width + j + 8)), dre_back8);
//		} // end-for
//		for (; j < width - 1; j += 1)
//		{
//			int com1 = smoothImg[(i + 1) * width + j + 1] - smoothImg[(i - 1) * width + j - 1];
//			int com2 = smoothImg[(i - 1) * width + j + 1] - smoothImg[(i + 1) * width + j - 1];
//			int gx;
//			int gy;
//			gx = abs(com1 + com2 + (smoothImg[i * width + j + 1] - smoothImg[i * width + j - 1]));
//			gy = abs(com1 - com2 + (smoothImg[(i + 1) * width + j] - smoothImg[(i - 1) * width + j]));
//			int sum;//梯度幅值
//			//sumFlag:梯度幅值计算方法
//			sum = gx + gy;
//			GradImg[i * width + j] = sum;
//			if (gx >= gy) DireImg[i * width + j] = 0;//垂直方向
//			else          DireImg[i * width + j] = 1;//水平方向
//		}
//	} // end-for
//#else
//	for (int i = 1; i < height - 1; i++) {
//		for (int j = 1; j < width - 1; j++) {
//			int com1 = smoothImg[(i + 1) * width + j + 1] - smoothImg[(i - 1) * width + j - 1];
//			int com2 = smoothImg[(i - 1) * width + j + 1] - smoothImg[(i + 1) * width + j - 1];
//			int gx;
//			int gy;
//			//根据不同方法计算梯度图
//			gx = abs(com1 + com2 + (smoothImg[i * width + j + 1] - smoothImg[i * width + j - 1]));
//			gy = abs(com1 - com2 + (smoothImg[(i + 1) * width + j] - smoothImg[(i - 1) * width + j]));
//			int sum;//梯度幅值
//			//sumFlag:梯度幅值计算方法
//			sum = gx + gy;
//			int index = i * width + j;//当前像素位置（data）
//			GradImg[index] = sum;
//			//梯度阈值与边缘方向
//			if (sum >= 1) {
//				if (gx >= gy) DireImg[index] = 0;//垂直方向
//				else          DireImg[index] = 1;//水平方向
//			} //end-if
//		} // end-for
//	} // end-for
//#endif
//}

//提取锚点
void ED::ComputeAnchorPoints()
{
	//memset(edgeImg, 0, width*height);
	for (int i = 2; i < height - 2; i++) {
		int start = 2;
		int inc = 1;
		//q取样间隔
		if (i % scanInterval != 0) { start = scanInterval; inc = scanInterval; }

		for (int j = start; j < width - 2; j += inc) {
			if (gradImg[i * width + j] < gradThresh) continue;//梯度值大于阈值继续
			//根据边缘方向计算梯度的差值
			if (dirImg[i * width + j] == EDGE_VERTICAL) {
				// 垂直方向，当前像素与上下像素的差值，两个差值都大于锚点阈值，则该像素为锚点
				int diff1 = gradImg[i * width + j] - gradImg[i * width + j - 1];
				int diff2 = gradImg[i * width + j] - gradImg[i * width + j + 1];
				if (diff1 >= anchorThresh && diff2 >= anchorThresh) {
					edgeImg[i * width + j] = ANCHOR_PIXEL;//edgeImg 确定锚点 ANCHOR_PIXEL = 254
					anchorPoints.push_back(Point(j, i));
				}
			}
			else {
				// 水平方向的像素
				int diff1 = gradImg[i * width + j] - gradImg[(i - 1) * width + j];
				int diff2 = gradImg[i * width + j] - gradImg[(i + 1) * width + j];
				if (diff1 >= anchorThresh && diff2 >= anchorThresh) {
					edgeImg[i * width + j] = ANCHOR_PIXEL;
					anchorPoints.push_back(Point(j, i));
				}
			} // end-else
		} //end-for-inner
	} //end-for-outer

	anchorNos = anchorPoints.size(); // 锚点总数
}


//连接锚点
void ED::JoinAnchorPointsUsingSortedAnchors()
{
	int* chainNos = new int[(width + height) * 8];//链接段数
	Point* pixels = new Point[width * height];//点
	StackNode* stack = new StackNode[width * height];//节点
	Chain* chains = new Chain[width * height];//链

	// A：梯度值递增，相同梯度值坐标次序递减，容量为锚点数量，A记录像素位置
	int* A = sortAnchorsByGradValue1();
	int totalPixels = 0; //总像素点（）

	//从梯度值最大的锚点开始连接
	for (int k = anchorNos - 1; k >= 0; k--) {
		int pixelOffset = A[k]; //起点像素位置

		int i = pixelOffset / width;//像素 行
		int j = pixelOffset % width;//像素 列

		if (edgeImg[i * width + j] != ANCHOR_PIXEL) continue;

		//初始化chains[0]
		chains[0].len = 0;
		chains[0].parent = -1;
		chains[0].dir = 0;
		chains[0].children[0] = chains[0].children[1] = -1;
		chains[0].pixels = NULL;

		int noChains = 1;//链编号
		int len = 0;//链中像素数量
		int duplicatePixelCount = 0;//重复像素点计数
		int top = -1;  // top of the stack 栈的顶

		//垂直方向的像素
		if (dirImg[i * width + j] == EDGE_VERTICAL) {
			//top = 0
			stack[++top].r = i;//行
			stack[top].c = j;//列
			stack[top].dir = DOWN;//方向向下连接
			stack[top].parent = 0;//父

			//top = 1
			stack[++top].r = i;
			stack[top].c = j;
			stack[top].dir = UP;//向上连接
			stack[top].parent = 0;
		}
		//水平
		else {
			//top = 0
			stack[++top].r = i;
			stack[top].c = j;
			stack[top].dir = RIGHT;//向右连接
			stack[top].parent = 0;

			//top = 1
			stack[++top].r = i;
			stack[top].c = j;
			stack[top].dir = LEFT;//向左连接
			stack[top].parent = 0;
		} //end-else

		  // While the stack is not empty：链不能在增长
	StartOfWhile: //
		//两个方向依次连接
		while (top >= 0) {
			int r = stack[top].r;//行
			int c = stack[top].c;//列
			int dir = stack[top].dir;
			int parent = stack[top].parent;
			top--;

			if (edgeImg[r * width + c] != EDGE_PIXEL) duplicatePixelCount++;//当前点不为255 则计数，记录重复的点数，一个像素点会有两个方向检索状态

			//初始化chains[noChains]
			chains[noChains].dir = dir;   // traversal direction
			chains[noChains].parent = parent;
			chains[noChains].children[0] = chains[noChains].children[1] = -1;
			chains[noChains].pixels = &pixels[len];  //len初始值为 0 ，链中像素的数量。最后一个像素坐标的引用

			int chainLen = 0; //链长

			//记录链中像素的坐标
			pixels[len].y = r;//起点
			pixels[len].x = c;
			len++;//链中像素的数量
			chainLen++;//链长

			if (dir == LEFT) {
				//当前像素是否是水平的，下一个点是水平则只能向左检索，右边是上一个点，不是水平则跳出循环
				while (dirImg[r * width + c] == EDGE_HORIZONTAL) {
					edgeImg[r * width + c] = EDGE_PIXEL; //当前像素=255
					//   A
					//   B x 
					//   C 
					//
					// 锚点上下像素置为0
					if (edgeImg[(r - 1) * width + c] == ANCHOR_PIXEL) edgeImg[(r - 1) * width + c] = 0;
					if (edgeImg[(r + 1) * width + c] == ANCHOR_PIXEL) edgeImg[(r + 1) * width + c] = 0;

					// 看邻域中是否有边缘像素：254、255// 寻找下一个点的坐标
					if (edgeImg[r * width + c - 1] >= ANCHOR_PIXEL) { c--; }//左边像素是否为锚点或边缘点，是则坐标移到该锚点或边缘点
					else if (edgeImg[(r - 1) * width + c - 1] >= ANCHOR_PIXEL) { r--; c--; }//左上角
					else if (edgeImg[(r + 1) * width + c - 1] >= ANCHOR_PIXEL) { r++; c--; }//左下角
					//左边三像素均不为锚点或边缘点，则找三点梯度的最大值
					else {
						int A = gradImg[(r - 1) * width + c - 1];
						int B = gradImg[r * width + c - 1];
						int C = gradImg[(r + 1) * width + c - 1];

						if (A > B) {
							if (A > C) r--;//A最大，坐标移到A
							else       r++;//C最大
						}
						else  if (C > B) r++;//C最大，坐标移到C
						c--;
					} //end-else

					//移到下一个点，处理新点，若这个新点 不 可以连接上一个点：
					//如果该点为边缘点或梯度值小于梯度阈值（梯度阈值筛选边缘候选点，锚点边缘点均在候选点中），则这个新点 不 可以连接上一个点：
					if (edgeImg[r * width + c] == EDGE_PIXEL || gradImg[r * width + c] < gradThresh) {
						if (chainLen > 0) {
							chains[noChains].len = chainLen;//chainLen:491行，chainLen++，当前编号的链的长度
							chains[parent].children[0] = noChains;//当前链的编号记录到父链的子链中，左边就是children[0]，右边1
							noChains++;//当前链结束，开始下一段链（检索到边缘以及换向都会导致链的结束）
						}
						goto StartOfWhile;//检索到边缘后，从另一边开始检索连接，469行，若两边都检索完了，则当前锚点的连接结束
					}

					//若这个新点可以连接上一个点：
					pixels[len].y = r;//记录新点的坐标，len：当前链的像素总数
					pixels[len].x = c;
					len++; //更新len，此时len比像素总数大1
					chainLen++;
				} //end-while，返回497行的while判断新点的方向，若为水平则继续循环，若不是，进行下一步

				//当新像素点不是水平方向：更新top，向垂直方向检索，且新点不能作为上一个链的点并计数
				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = DOWN;
				stack[top].parent = noChains;

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = UP;
				stack[top].parent = noChains;

				len--;
				chainLen--;

				chains[noChains].len = chainLen;//当前链的像素数=长度
				chains[parent].children[0] = noChains; //当前链的编号记录到父链的子链中，左边就是children[0]，右边1
				noChains++;//更新链的编号，当前链结束，开始下一段链（检索到边缘以及换向都会导致链的结束）

			}
			else if (dir == RIGHT) {
				while (dirImg[r * width + c] == EDGE_HORIZONTAL) {
					edgeImg[r * width + c] = EDGE_PIXEL;

					// The edge is horizontal. Look RIGHT
					//
					//     A
					//   x B
					//     C
					//
					// cleanup up&down pixels
					if (edgeImg[(r + 1) * width + c] == ANCHOR_PIXEL) edgeImg[(r + 1) * width + c] = 0;
					if (edgeImg[(r - 1) * width + c] == ANCHOR_PIXEL) edgeImg[(r - 1) * width + c] = 0;

					// Look if there is an edge pixel in the neighbors
					if (edgeImg[r * width + c + 1] >= ANCHOR_PIXEL) { c++; }
					else if (edgeImg[(r + 1) * width + c + 1] >= ANCHOR_PIXEL) { r++; c++; }
					else if (edgeImg[(r - 1) * width + c + 1] >= ANCHOR_PIXEL) { r--; c++; }
					else {
						// else -- follow max. pixel to the RIGHT
						int A = gradImg[(r - 1) * width + c + 1];
						int B = gradImg[r * width + c + 1];
						int C = gradImg[(r + 1) * width + c + 1];

						if (A > B) {
							if (A > C) r--;       // A
							else       r++;       // C
						}
						else if (C > B) r++;  // C
						c++;
					} //end-else

					if (edgeImg[r * width + c] == EDGE_PIXEL || gradImg[r * width + c] < gradThresh) {
						if (chainLen > 0) {
							chains[noChains].len = chainLen;
							chains[parent].children[1] = noChains;
							noChains++;
						} // end-if
						goto StartOfWhile;
					} //end-else


					pixels[len].y = r;
					pixels[len].x = c;
					len++;
					chainLen++;
				} //end-while

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = DOWN;  // Go down
				stack[top].parent = noChains;

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = UP;   // Go up
				stack[top].parent = noChains;

				len--;
				chainLen--;

				chains[noChains].len = chainLen;
				chains[parent].children[1] = noChains;
				noChains++;

			}
			else if (dir == UP) {
				while (dirImg[r * width + c] == EDGE_VERTICAL) {
					edgeImg[r * width + c] = EDGE_PIXEL;

					// The edge is vertical. Look UP
					//
					//   A B C
					//     x
					//
					// Cleanup left & right pixels
					if (edgeImg[r * width + c - 1] == ANCHOR_PIXEL) edgeImg[r * width + c - 1] = 0;
					if (edgeImg[r * width + c + 1] == ANCHOR_PIXEL) edgeImg[r * width + c + 1] = 0;

					// Look if there is an edge pixel in the neighbors
					if (edgeImg[(r - 1) * width + c] >= ANCHOR_PIXEL) { r--; }
					else if (edgeImg[(r - 1) * width + c - 1] >= ANCHOR_PIXEL) { r--; c--; }
					else if (edgeImg[(r - 1) * width + c + 1] >= ANCHOR_PIXEL) { r--; c++; }
					else {
						// else -- follow the max. pixel UP
						int A = gradImg[(r - 1) * width + c - 1];
						int B = gradImg[(r - 1) * width + c];
						int C = gradImg[(r - 1) * width + c + 1];

						if (A > B) {
							if (A > C) c--;
							else       c++;
						}
						else if (C > B) c++;
						r--;
					} //end-else

					if (edgeImg[r * width + c] == EDGE_PIXEL || gradImg[r * width + c] < gradThresh) {
						if (chainLen > 0) {
							chains[noChains].len = chainLen;
							chains[parent].children[0] = noChains;
							noChains++;
						} // end-if
						goto StartOfWhile;
					} //end-else


					pixels[len].y = r;
					pixels[len].x = c;

					len++;
					chainLen++;
				} //end-while

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = RIGHT;
				stack[top].parent = noChains;

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = LEFT;
				stack[top].parent = noChains;

				len--;
				chainLen--;

				chains[noChains].len = chainLen;
				chains[parent].children[0] = noChains;
				noChains++;

			}
			else { // dir == DOWN
				while (dirImg[r * width + c] == EDGE_VERTICAL) {
					edgeImg[r * width + c] = EDGE_PIXEL;

					// The edge is vertical
					//
					//     x
					//   A B C
					//
					// cleanup side pixels
					if (edgeImg[r * width + c + 1] == ANCHOR_PIXEL) edgeImg[r * width + c + 1] = 0;
					if (edgeImg[r * width + c - 1] == ANCHOR_PIXEL) edgeImg[r * width + c - 1] = 0;

					// Look if there is an edge pixel in the neighbors
					if (edgeImg[(r + 1) * width + c] >= ANCHOR_PIXEL) { r++; }
					else if (edgeImg[(r + 1) * width + c + 1] >= ANCHOR_PIXEL) { r++; c++; }
					else if (edgeImg[(r + 1) * width + c - 1] >= ANCHOR_PIXEL) { r++; c--; }
					else {
						// else -- follow the max. pixel DOWN
						int A = gradImg[(r + 1) * width + c - 1];
						int B = gradImg[(r + 1) * width + c];
						int C = gradImg[(r + 1) * width + c + 1];

						if (A > B) {
							if (A > C) c--;       // A
							else       c++;       // C
						}
						else if (C > B) c++;  // C
						r++;
					} //end-else

					if (edgeImg[r * width + c] == EDGE_PIXEL || gradImg[r * width + c] < gradThresh) {
						if (chainLen > 0) {
							chains[noChains].len = chainLen;
							chains[parent].children[1] = noChains;
							noChains++;
						} // end-if
						goto StartOfWhile;
					} //end-else

					pixels[len].y = r;
					pixels[len].x = c;

					len++;
					chainLen++;
				} //end-while

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = RIGHT;
				stack[top].parent = noChains;

				stack[++top].r = r;
				stack[top].c = c;
				stack[top].dir = LEFT;
				stack[top].parent = noChains;

				len--;
				chainLen--;

				chains[noChains].len = chainLen;
				chains[parent].children[1] = noChains;
				noChains++;
			} //end-else		
		} //end-while//一直检索到停止条件：非边缘候选点、以及以及标记的边缘点时，当前锚点的链结束（两边都达到停止条件）

		//去除短的链，阈值：minPathLen
		if (len - duplicatePixelCount < minPathLen) {
			for (int k = 0; k < len; k++) {

				edgeImg[pixels[k].y * width + pixels[k].x] = 0;
				edgeImg[pixels[k].y * width + pixels[k].x] = 0;
			} //end-for
		}
		//长度符合的链
		else {

			int noSegmentPixels = 0;
			//返回长度最长的链
			int totalLen = LongestChain(chains, chains[0].children[1]);

			if (totalLen > 0) {
				// Retrieve the chainNos
				int count = RetrieveChainNos(chains, chains[0].children[1], chainNos);

				// Copy these pixels in the reverse order
				for (int k1 = count - 1; k1 >= 0; k1--) {
					int chainNo = chainNos[k1];

#if 1
					/* See if we can erase some pixels from the last chain. This is for cleanup */

					int fr = chains[chainNo].pixels[chains[chainNo].len - 1].y;
					int fc = chains[chainNo].pixels[chains[chainNo].len - 1].x;

					int index = noSegmentPixels - 2;
					while (index >= 0) {
						int dr = abs(fr - segmentPoints[segmentNos][index].y);
						int dc = abs(fc - segmentPoints[segmentNos][index].x);

						if (dr <= 1 && dc <= 1) {
							// neighbors. Erase last pixel
							segmentPoints[segmentNos].pop_back();
							noSegmentPixels--;
							index--;
						}
						else break;
					} //end-while

					if (chains[chainNo].len > 1 && noSegmentPixels > 0) {
						fr = chains[chainNo].pixels[chains[chainNo].len - 2].y;
						fc = chains[chainNo].pixels[chains[chainNo].len - 2].x;

						int dr = abs(fr - segmentPoints[segmentNos][noSegmentPixels - 1].y);
						int dc = abs(fc - segmentPoints[segmentNos][noSegmentPixels - 1].x);

						if (dr <= 1 && dc <= 1) chains[chainNo].len--;
					} //end-if
#endif

					for (int l = chains[chainNo].len - 1; l >= 0; l--) {
						segmentPoints[segmentNos].push_back(chains[chainNo].pixels[l]);
						noSegmentPixels++;
					} //end-for

					chains[chainNo].len = 0;  // Mark as copied
				} //end-for
			} //end-if

			totalLen = LongestChain(chains, chains[0].children[0]);
			if (totalLen > 1) {
				// Retrieve the chainNos
				int count = RetrieveChainNos(chains, chains[0].children[0], chainNos);

				// Copy these chains in the forward direction. Skip the first pixel of the first chain
				// due to repetition with the last pixel of the previous chain
				int lastChainNo = chainNos[0];
				chains[lastChainNo].pixels++;
				chains[lastChainNo].len--;

				for (int k2 = 0; k2 < count; k2++) {
					int chainNo = chainNos[k2];

#if 1
					/* See if we can erase some pixels from the last chain. This is for cleanup */
					int fr = chains[chainNo].pixels[0].y;
					int fc = chains[chainNo].pixels[0].x;

					int index = noSegmentPixels - 2;
					while (index >= 0) {
						int dr = abs(fr - segmentPoints[segmentNos][index].y);
						int dc = abs(fc - segmentPoints[segmentNos][index].x);

						if (dr <= 1 && dc <= 1) {
							// neighbors. Erase last pixel
							segmentPoints[segmentNos].pop_back();
							noSegmentPixels--;
							index--;
						}
						else break;
					} //end-while

					int startIndex = 0;
					int chainLen = chains[chainNo].len;
					if (chainLen > 1 && noSegmentPixels > 0) {
						int fr = chains[chainNo].pixels[1].y;
						int fc = chains[chainNo].pixels[1].x;

						int dr = abs(fr - segmentPoints[segmentNos][noSegmentPixels - 1].y);
						int dc = abs(fc - segmentPoints[segmentNos][noSegmentPixels - 1].x);

						if (dr <= 1 && dc <= 1) { startIndex = 1; }
					} //end-if
#endif

					  /* Start a new chain & copy pixels from the new chain */
					for (int l = startIndex; l < chains[chainNo].len; l++) {
						segmentPoints[segmentNos].push_back(chains[chainNo].pixels[l]);
						noSegmentPixels++;
					} //end-for

					chains[chainNo].len = 0;  // Mark as copied
				} //end-for
			} //end-if


			  // See if the first pixel can be cleaned up
			int fr = segmentPoints[segmentNos][1].y;
			int fc = segmentPoints[segmentNos][1].x;


			int dr = abs(fr - segmentPoints[segmentNos][noSegmentPixels - 1].y);
			int dc = abs(fc - segmentPoints[segmentNos][noSegmentPixels - 1].x);


			if (dr <= 1 && dc <= 1) {
				segmentPoints[segmentNos].erase(segmentPoints[segmentNos].begin());
				noSegmentPixels--;
			} //end-if

			segmentNos++;
			segmentPoints.push_back(vector<Point>()); // create empty vector of points for segments

													  // Copy the rest of the long chains here
			for (int k3 = 2; k3 < noChains; k3++) {
				if (chains[k3].len < 2) continue;

				totalLen = LongestChain(chains, k3);

				if (totalLen >= 10) {

					// Retrieve the chainNos
					int count = RetrieveChainNos(chains, k3, chainNos);

					// Copy the pixels
					noSegmentPixels = 0;
					for (int k4 = 0; k4 < count; k4++) {
						int chainNo = chainNos[k4];

#if 1					
						/* See if we can erase some pixels from the last chain. This is for cleanup */
						int fr = chains[chainNo].pixels[0].y;
						int fc = chains[chainNo].pixels[0].x;

						int index = noSegmentPixels - 2;
						while (index >= 0) {
							int dr = abs(fr - segmentPoints[segmentNos][index].y);
							int dc = abs(fc - segmentPoints[segmentNos][index].x);

							if (dr <= 1 && dc <= 1) {
								// neighbors. Erase last pixel
								segmentPoints[segmentNos].pop_back();
								noSegmentPixels--;
								index--;
							}
							else break;
						} //end-while

						int startIndex = 0;
						int chainLen = chains[chainNo].len;
						if (chainLen > 1 && noSegmentPixels > 0) {
							int fr = chains[chainNo].pixels[1].y;
							int fc = chains[chainNo].pixels[1].x;

							int dr = abs(fr - segmentPoints[segmentNos][noSegmentPixels - 1].y);
							int dc = abs(fc - segmentPoints[segmentNos][noSegmentPixels - 1].x);

							if (dr <= 1 && dc <= 1) { startIndex = 1; }
						} //end-if
#endif
						  /* Start a new chain & copy pixels from the new chain */
						for (int l = startIndex; l < chains[chainNo].len; l++) {
							segmentPoints[segmentNos].push_back(chains[chainNo].pixels[l]);
							noSegmentPixels++;
						} //end-for

						chains[chainNo].len = 0;  // Mark as copied
					} //end-for
					segmentPoints.push_back(vector<Point>()); // create empty vector of points for segments
					segmentNos++;
				} //end-if          
			} //end-for

		} //end-else

	} //end-for-outer

	// pop back last segment from vector
	// because of one preallocation in the beginning, it will always empty
	segmentPoints.pop_back();

	// Clean up
	delete[] A;
	delete[] chains;
	delete[] stack;
	delete[] chainNos;
	delete[] pixels;
}

void ED::sortAnchorsByGradValue()
{
	auto sortFunc = [&](const Point& a, const Point& b)
	{
		return gradImg[a.y * width + a.x] > gradImg[b.y * width + b.x];
	};

	std::sort(anchorPoints.begin(), anchorPoints.end(), sortFunc);

	/*
	ofstream myFile;
	myFile.open("anchorsNew.txt");
	for (int i = 0; i < anchorPoints.size(); i++) {
		int x = anchorPoints[i].x;
		int y = anchorPoints[i].y;

		myFile << i << ". value: " << gradImg[y*width + x] << "  Cord: (" << x << "," << y << ")" << endl;
	}
	myFile.close();


	vector<Point> temp(anchorNos);

	int x, y, i = 0;
	char c;
	std::ifstream infile("cords.txt");
	while (infile >> x >> c >> y && c == ',') {
		temp[i] = Point(x, y);
		i++;
	}

	anchorPoints = temp;
	*/
}

//根据梯度值排序锚点，A：梯度值递增，相同梯度值坐标次序递减，容量为锚点数量
int* ED::sortAnchorsByGradValue1()
{
	int SIZE = 128 * 256; //原参数：128*256
	int* C = new int[SIZE];
	memset(C, 0, sizeof(int) * SIZE);//int-4字节，初始化每个字节为0（00000000）并拷贝到C中

	// Count the number of grad values统计锚点处 梯度值的数量
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {
			if (edgeImg[i * width + j] != ANCHOR_PIXEL) continue;

			int grad = gradImg[i * width + j];//锚点处的梯度值
			C[grad]++;//该梯度值所在位置计数
		} //end-for
	} //end-for 

	// Compute indices指数，递增
	for (int i = 1; i < SIZE; i++) C[i] += C[i - 1];

	int noAnchors = C[SIZE - 1];//梯度值的数量
	int* A = new int[noAnchors]; //A：锚点数量？
	memset(A, 0, sizeof(int) * noAnchors);


	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {
			if (edgeImg[i * width + j] != ANCHOR_PIXEL) continue;

			int grad = gradImg[i * width + j];
			int index = --C[grad]; //index：对应梯度值指数-1
			A[index] = i * width + j;    // anchor's offset 锚点的data坐标
		} //end-for
	} //end-for  
	//cout << "11111111111   111   " << C[SIZE - 1] << endl;
	delete[] C;


	//ofstream myFile;
	//myFile.open("aNew.txt");
	//for (int i = 0; i < noAnchors; i++)
	//	myFile << A[i] << endl;
	//myFile.close(); 

	//A：梯度值递增，相同梯度值坐标次序递减，容量为锚点数量
	return A;
}

//返回路径最长的链
int ED::LongestChain(Chain* chains, int root) {
	if (root == -1 || chains[root].len == 0) return 0;

	int len0 = 0;
	if (chains[root].children[0] != -1) len0 = LongestChain(chains, chains[root].children[0]);//

	int len1 = 0;
	if (chains[root].children[1] != -1) len1 = LongestChain(chains, chains[root].children[1]);

	int max = 0;

	if (len0 >= len1) {
		max = len0;
		chains[root].children[1] = -1;

	}
	else {
		max = len1;
		chains[root].children[0] = -1;
	} //end-else

	return chains[root].len + max;
} //end-LongestChain
//链计数并返回对应的编号
int ED::RetrieveChainNos(Chain* chains, int root, int chainNos[])
{
	int count = 0;

	while (root != -1) {
		chainNos[count] = root;
		count++;

		if (chains[root].children[0] != -1) root = chains[root].children[0];
		else                                root = chains[root].children[1];
	} //end-while

	return count;
}
