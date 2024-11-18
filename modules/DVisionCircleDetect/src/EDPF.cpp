#include "EDPF.h"

using namespace cv;
using namespace std;

EDPF::EDPF(Mat srcImage)
	:ED(srcImage, PREWITT_OPERATOR, 11, 3)//初始化，
{
	// Validate Edge Segments
	//sigma /= 2.5;
	//GaussianBlur(srcImage, smoothImage, Size(), sigma, sigma); // calculate kernel from sigma

	validateEdgeSegments();
}

void EDPF::validateEdgeSegments()
{
	divForTestSegment = 2.25; // Some magic number :-)2.25
	memset(edgeImg, 0, width * height); // clear edge image

	H = new double[MAX_GRAD_VALUE];//梯度范数累计分布函数
	memset(H, 0, sizeof(double) * MAX_GRAD_VALUE);

	//返回图像梯度图，计算梯度概率函数（递减）
	gradImg = ComputePrewitt3x3();
	np = 0;
	for (int i = 0; i < segmentNos; i++) {
		int len = segmentPoints[i].size();
		np += (len * (len - 1)) / 2;
	} //end-for

	// Validate segments
	for (int i = 0; i < segmentNos; i++) {
		TestSegment(i, 0, segmentPoints[i].size() - 1);
	} //end-for

	ExtractNewSegments();

	// clean space		  
	delete[] H;
	delete[] gradImg;
}

//返回图像梯度图，计算梯度概率函数（递减）
inline
short* EDPF::ComputePrewitt3x3()
{
	short* gradImg = new short[width * height];
	memset(gradImg, 0, sizeof(short) * width * height);

	int* grads = new int[MAX_GRAD_VALUE];
	memset(grads, 0, sizeof(int) * MAX_GRAD_VALUE);

	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width - 1; j++) {
		
			int com1 = srcImg[(i + 1) * width + j + 1] - srcImg[(i - 1) * width + j - 1];
			int com2 = srcImg[(i - 1) * width + j + 1] - srcImg[(i + 1) * width + j - 1];

			int gx = abs(com1 + com2 + (srcImg[i * width + j + 1] - srcImg[i * width + j - 1]));
			int gy = abs(com1 - com2 + (srcImg[(i + 1) * width + j] - srcImg[(i - 1) * width + j]));

			int g = gx + gy;

			gradImg[i * width + j] = g;
			grads[g]++;
		} // end-for
	} //end-for

	 // Compute probability function H
	int size = (width - 2) * (height - 2);

	for (int i = MAX_GRAD_VALUE - 1; i > 0; i--)
		grads[i - 1] += grads[i];

	for (int i = 0; i < MAX_GRAD_VALUE; i++)
		H[i] = (double)grads[i] / ((double)size);

	delete[] grads;
	return gradImg;
}


void EDPF::TestSegment(int i, int index1, int index2)
{

	int chainLen = index2 - index1 + 1;
	if (chainLen < minPathLen)
		return;

	int minGrad = 1 << 30;//2^30
	int minGradIndex;
	for (int k = index1; k <= index2; k++) {
		int r = segmentPoints[i][k].y;
		int c = segmentPoints[i][k].x;
		if (gradImg[r * width + c] < minGrad) { minGrad = gradImg[r * width + c]; minGradIndex = k; }//第一个最小值
	} //end-for

	 // Compute nfa,计算 number of false alarms
	double nfa = NFA(H[minGrad], (int)(chainLen / divForTestSegment));//divForTestSegment = 2.25

	if (nfa <= EPSILON) {
		for (int k = index1; k <= index2; k++) {
			int r = segmentPoints[i][k].y;
			int c = segmentPoints[i][k].x;

			edgeImg[r * width + c] = 255;
		} //end-for

		return;
	} //end-if  

	// Split into two halves. We divide at the point where the gradient is the minimum
	int end = minGradIndex - 1;
	while (end > index1) {
		int r = segmentPoints[i][end].y;
		int c = segmentPoints[i][end].x;

		if (gradImg[r * width + c] <= minGrad) end--;
		else break;
	} //end-while

	int start = minGradIndex + 1;
	while (start < index2) {
		int r = segmentPoints[i][start].y;
		int c = segmentPoints[i][start].x;

		if (gradImg[r * width + c] <= minGrad) start++;
		else break;
	} //end-while

	TestSegment(i, index1, end);
	TestSegment(i, start, index2);
}


void EDPF::ExtractNewSegments()
{
	//vector<Point> *segments = &segmentPoints[segmentNos];
	vector< vector<Point> > validSegments;
	int noSegments = 0;

	for (int i = 0; i < segmentNos; i++) {
		int start = 0;
		while (start < segmentPoints[i].size()) {

			while (start < segmentPoints[i].size()) {
				int r = segmentPoints[i][start].y;
				int c = segmentPoints[i][start].x;

				if (edgeImg[r * width + c]) break;
				start++;
			} //end-while

			int end = start + 1;
			while (end < segmentPoints[i].size()) {
				int r = segmentPoints[i][end].y;
				int c = segmentPoints[i][end].x;

				if (edgeImg[r * width + c] == 0) break;
				end++;
			} //end-while

			int len = end - start;
			if (len >= 10) {
				validSegments.push_back(vector<Point>());
				vector<Point> subVec(&segmentPoints[i][start], &segmentPoints[i][end - 1]);
				validSegments[noSegments] = subVec;
				noSegments++;
			} //end-else

			start = end + 1;
		} //end-while
	} //end-for

	 // Copy to ed
	segmentPoints = validSegments;

	segmentNos = noSegments;
}

//---------------------------------------------------------------------------
// Number of false alarms code as suggested by Desolneux, Moisan and Morel (DMM)
//
double EDPF::NFA(double prob, int len)
{
	double nfa = np;
	for (int i = 0; i<len && nfa > EPSILON; i++)
		nfa *= prob;

	return nfa;
}