#pragma once

#include "ED.h"

#define MAX_GRAD_VALUE 6*256  //最大梯度值
#define EPSILON 1.0 //ε

class EDPF : public ED {
public:
	EDPF(cv::Mat srcImage);
	EDPF() {};

private:
	double divForTestSegment;
	double* H;
	int np;
	short* gradImg;

	void validateEdgeSegments();
	short* ComputePrewitt3x3(); // differs from base class's prewit function (calculates H)
	void TestSegment(int i, int index1, int index2);
	void ExtractNewSegments();
	double NFA(double prob, int len);
};


