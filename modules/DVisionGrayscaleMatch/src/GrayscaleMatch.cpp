#include "GrayscaleMatch.h"
#include <opencv2/core/simd_intrinsics.hpp>
#include<thread>


#define VISION_TOLERANCE 0.0000001
#define D2R (CV_PI / 180.0)
#define R2D (180.0 / CV_PI)
#define MATCH_CANDIDATE_NUM 5

using namespace std;
using namespace cv;


void GrayscaleMatch::setTempimage(cv::Mat _temp, bool useMask, cv::Mat _mask)
{
	m_params.temp = _temp;
	m_params.useMask = useMask;
	if(useMask)
		m_params.mask = _mask;
}

void GrayscaleMatch::setScore(float score)
{
	m_params.score = score;
}

void GrayscaleMatch::setMaxtargs(float maxtargs)
{
	m_params.maxtargs = maxtargs;
}

void GrayscaleMatch::setPyramidLayer(float pyramidLayer)
{
	m_params.pyramidLayer = pyramidLayer;
}

void GrayscaleMatch::setAngle(float angleStart, float angleRange)
{
	m_params.angleStart = angleStart;
	m_params.angleRange = angleRange;
}

void GrayscaleMatch::setMaxOverlap(float maxOverlap)
{
	m_params.maxOverlap = maxOverlap;
}

void GrayscaleMatch::setPolarity(bool polarity)
{
	m_params.polarity = polarity;
}

void GrayscaleMatch::setSubpixel(bool subpixel)
{
	m_params.subpixel = subpixel;
}

void GrayscaleMatch::setMean(bool useMean, float mean)
{
	m_params.useMean = useMean;
	m_params.mean = mean;
}

void GrayscaleMatch::setSDV(bool useSDV, float SDV)
{
	m_params.useSDV = useSDV;
	m_params.SDV = SDV;
}


// 初始化
bool GrayscaleMatch::initialize()
{
	// 模板图
	if (m_params.temp.empty())
	{
		std::cout << "template image is empty!" << endl;
		return false;
	}
	switch (m_params.temp.channels())
	{
	case 3:
		cvtColor(m_params.temp, m_params.temp, cv::COLOR_BGR2GRAY);
		break;
	case 4:
		cvtColor(m_params.temp, m_params.temp, cv::COLOR_BGRA2GRAY);
		break;
	default:
		break;
	}
	cout << m_params.temp.size() << endl;
	// 掩码图
	if (m_params.useMask)
	{
		if (m_params.mask.empty())
		{
			std::cout << "mask image is empty!" << endl;
			return false;
		}
		if (m_params.mask.size() != m_params.temp.size())
		{
			std::cout << "the size of mask image is not equal to the template image!" << endl;
			return false; 
		}
		switch (m_params.mask.channels())
		{
		case 3:
			cvtColor(m_params.mask, m_params.mask, cv::COLOR_BGR2GRAY);
			break;
		case 4:
			cvtColor(m_params.mask, m_params.mask, cv::COLOR_BGRA2GRAY);
			break;
		default:
			break;
		}
		threshold(m_params.mask, m_params.mask, 0/*threshold*/, 1.0/*maxVal*/, THRESH_BINARY);
	}

	// 灰度图
	if (m_gray.empty())
	{
		std::cout << "src image is empty!" << endl;
		return false;
	}
	switch (m_gray.channels())
	{
	case 3:
		cvtColor(m_gray, m_gray, cv::COLOR_BGR2GRAY);
		break;
	case 4:
		cvtColor(m_gray, m_gray, cv::COLOR_BGRA2GRAY);
		break;
	default:
		break;
	}

	// 金字塔参数
	if (m_params.pyramidLayer >= 10)
	{
		ushort TopLayer = 0;
		ushort iMinDstLength = (int)cv::sqrt(m_params.pyramidLayer);
		ushort iMinReduceArea = iMinDstLength * iMinDstLength;
		int iArea = m_params.temp.cols * m_params.temp.rows;
		while (iArea > iMinReduceArea)
		{
			iArea /= 4;
			TopLayer++;
		}
		m_params.pyramidLayer = TopLayer;
	}
	int iTopLayer = m_params.pyramidLayer;

	// 旋转角度参数
	float len = 0.5 * sqrt(m_params.temp.rows * m_params.temp.rows + m_params.temp.cols * m_params.temp.cols);
	m_angleStep = 229.18 / len;
	m_angleStep = m_angleStep < FLT_MIN ? 1 : m_angleStep;
	m_anglePyr.resize(iTopLayer + 1);
	m_anglePyr[0] = m_angleStep;
	for (int i = 1; i < m_anglePyr.size(); i++)
	{
		m_anglePyr[i] = 2.0f * m_anglePyr[i - 1];
	}

	// 匹配分数
	m_scorePyr.resize(iTopLayer + 1);
	m_scorePyr[0] = m_params.score;
	for (int iLayer = 1; iLayer <= iTopLayer; iLayer++)
	{
		m_scorePyr[iLayer] = m_scorePyr[iLayer - 1] * 0.95;
	}

	m_TemplData.clear();
	if (m_params.useMask)
	{
		buildPyramid(m_params.mask, m_maskPyr, iTopLayer);
		for (int i = 0; i < m_maskPyr.size(); i++)
		{
			threshold(m_maskPyr[i], m_maskPyr[i], 0/*threshold*/, 1.0/*maxVal*/, THRESH_BINARY);
		}
	}

	buildPyramid(m_gray, m_grayPyr, iTopLayer);
	buildPyramid(m_params.temp, m_TemplData.vecPyramid, iTopLayer);
	m_TemplData.iBorderColor = mean(m_params.temp).val[0] < 128 ? 255 : 0;
	int iSize = m_TemplData.vecPyramid.size();
	m_TemplData.resize(iSize);
	for (int i = 0; i < iSize; i++)
	{
		double invArea = 1. / ((double)m_TemplData.vecPyramid[i].rows * m_TemplData.vecPyramid[i].cols);
		Scalar templMean, templSdv;
		double templNorm = 0, templSum2 = 0;

		meanStdDev(m_TemplData.vecPyramid[i], templMean, templSdv);
		templNorm = templSdv[0] * templSdv[0] + templSdv[1] * templSdv[1] + templSdv[2] * templSdv[2] + templSdv[3] * templSdv[3];

		if (templNorm < DBL_EPSILON)
		{
			m_TemplData.vecResultEqual1[i] = true;
		}
		templSum2 = templNorm + templMean[0] * templMean[0] + templMean[1] * templMean[1] + templMean[2] * templMean[2] + templMean[3] * templMean[3];


		templSum2 /= invArea;
		templNorm = std::sqrt(templNorm);
		templNorm /= std::sqrt(invArea); // care of accuracy here
		m_TemplData.vecInvArea[i] = invArea;
		m_TemplData.vecTemplMean[i] = templMean;
		m_TemplData.vecTemplNorm[i] = templNorm;
	}
	m_TemplData.bIsPatternLearned = true;

	return true;
}


inline int IM_Conv_CVSIMD(unsigned char* pCharKernel, unsigned char* pCharConv, int iLength)
{
	const int iBlockSize = CV_SIMD_WIDTH / sizeof(uchar);
	const int Block = iLength / iBlockSize;
	v_uint32 SumV = vx_setzero_u32();
	v_uint8 Zero = vx_setzero_u8();
	for (int Y = 0; Y < Block * iBlockSize; Y += iBlockSize)
	{
		v_uint8 SrcK = vx_load(pCharKernel + Y);
		v_uint8 SrcC = vx_load(pCharConv + Y);
		v_uint8 SrcK_L, SrcK_H, SrcC_L, SrcC_H;
		v_zip(SrcK, Zero, SrcK_L, SrcK_H);
		v_zip(SrcC, Zero, SrcC_L, SrcC_H);
		v_uint32 SumT1 = v_dotprod_expand(SrcK_L, SrcC_L);
		v_uint32 SumT2 = v_dotprod_expand(SrcK_H, SrcC_H);
		v_uint32 SumT = SumT1 + SumT2;
		SumV += SumT;
	}
	int Sum = v_reduce_sum(SumV);
	for (int Y = Block * iBlockSize; Y < iLength; Y++)
	{
		Sum += pCharKernel[Y] * pCharConv[Y];
	}
	return Sum;
}
inline int IM_Conv32f_CVSIMD(float* pCharKernel, float* pCharConv, int iLength)
{
	const int iBlockSize = CV_SIMD_WIDTH / sizeof(float);
	const int Block = iLength / iBlockSize;
	v_float32 SumV = vx_setzero_f32();
	for (int Y = 0; Y < Block * iBlockSize; Y += iBlockSize)
	{
		v_float32 SrcK = vx_load(pCharKernel + Y);
		v_float32 SrcC = vx_load(pCharConv + Y);
		v_float32 SumT = SrcK * SrcC;
		SumV += SumT;
	}
	int Sum = v_reduce_sum(SumV);
	for (int Y = Block * iBlockSize; Y < iLength; Y++)
	{
		Sum += pCharKernel[Y] * pCharConv[Y];
	}
	return Sum;
}

inline Point2f ptRotatePt2f(Point2f ptInput, Point2f ptOrg, double dAngle)
{
	double dWidth = ptOrg.x * 2;
	double dHeight = ptOrg.y * 2;
	double dY1 = dHeight - ptInput.y, dY2 = dHeight - ptOrg.y;

	double dX = (ptInput.x - ptOrg.x) * cos(dAngle) - (dY1 - ptOrg.y) * sin(dAngle) + ptOrg.x;
	double dY = (ptInput.x - ptOrg.x) * sin(dAngle) + (dY1 - ptOrg.y) * cos(dAngle) + dY2;

	dY = -dY + dHeight;
	return Point2f((float)dX, (float)dY);
}
cv::Size GetBestRotationSize(cv::Size sizeSrc, cv::Size sizeDst, double dRAngle)
{
	double dRAngle_radian = dRAngle * D2R;
	Point ptLT(0, 0), ptLB(0, sizeSrc.height - 1), ptRB(sizeSrc.width - 1, sizeSrc.height - 1), ptRT(sizeSrc.width - 1, 0);
	Point2f ptCenter((sizeSrc.width - 1) / 2.0f, (sizeSrc.height - 1) / 2.0f);
	Point2f ptLT_R = ptRotatePt2f(Point2f(ptLT), ptCenter, dRAngle_radian);
	Point2f ptLB_R = ptRotatePt2f(Point2f(ptLB), ptCenter, dRAngle_radian);
	Point2f ptRB_R = ptRotatePt2f(Point2f(ptRB), ptCenter, dRAngle_radian);
	Point2f ptRT_R = ptRotatePt2f(Point2f(ptRT), ptCenter, dRAngle_radian);

	float fTopY = max(max(ptLT_R.y, ptLB_R.y), max(ptRB_R.y, ptRT_R.y));
	float fBottomY = min(min(ptLT_R.y, ptLB_R.y), min(ptRB_R.y, ptRT_R.y));
	float fRightX = max(max(ptLT_R.x, ptLB_R.x), max(ptRB_R.x, ptRT_R.x));
	float fLeftX = min(min(ptLT_R.x, ptLB_R.x), min(ptRB_R.x, ptRT_R.x));


	if (dRAngle > 360)
		dRAngle -= 360;
	else if (dRAngle < 0)
		dRAngle += 360;

	if (fabs(fabs(dRAngle) - 90) < VISION_TOLERANCE || fabs(fabs(dRAngle) - 270) < VISION_TOLERANCE)
	{
		return Size(sizeSrc.height, sizeSrc.width);
	}
	else if (fabs(dRAngle) < VISION_TOLERANCE || fabs(fabs(dRAngle) - 180) < VISION_TOLERANCE)
	{
		return sizeSrc;
	}

	double dAngle = dRAngle;

	if (dAngle > 0 && dAngle < 90)
	{
		;
	}
	else if (dAngle > 90 && dAngle < 180)
	{
		dAngle -= 90;
	}
	else if (dAngle > 180 && dAngle < 270)
	{
		dAngle -= 180;
	}
	else if (dAngle > 270 && dAngle < 360)
	{
		dAngle -= 270;
	}

	float fH1 = sizeDst.width * sin(dAngle * D2R) * cos(dAngle * D2R);
	float fH2 = sizeDst.height * sin(dAngle * D2R) * cos(dAngle * D2R);

	int iHalfHeight = (int)ceil(fTopY - ptCenter.y - fH1);
	int iHalfWidth = (int)ceil(fRightX - ptCenter.x - fH2);
	Size sizeRet(iHalfWidth * 2, iHalfHeight * 2);

	bool bWrongSize = (sizeDst.width < sizeRet.width&& sizeDst.height > sizeRet.height)
			|| (sizeDst.width > sizeRet.width && sizeDst.height < sizeRet.height
			|| sizeDst.area() > sizeRet.area());
	if (bWrongSize)
		sizeRet = Size(int(fRightX - fLeftX + 0.5), int(fTopY - fBottomY + 0.5));
	return sizeRet;
}

void GrayscaleMatch::MatchTemplateMask(cv::Mat& matSrc, s_TemplData* pTemplData, cv::Mat& matResult, int iLayer)
{
	cv::Mat Templ = pTemplData->vecPyramid[iLayer];
	cv::Mat Mask = m_maskPyr[iLayer];
	Templ.convertTo(Templ, CV_32F);
	Mask.convertTo(Mask, CV_32F);
	Size corrSize(matSrc.cols - Templ.cols + 1, matSrc.rows - Templ.rows + 1);
	Mat Img; matSrc.convertTo(Img, CV_32F);
	matResult.create(corrSize, CV_32F);
	matResult.setTo(0);
	Mat img_mask = Mat(corrSize, CV_32FC1, Scalar(0));
	Scalar mask_sum = sum(Mask);
	Mat templx_mask = Mask.mul(Mask.mul(Templ - sum(Mask.mul(Templ)).div(mask_sum)));


	int  t_r_end = Templ.rows, t_r = 0;
	for (int r = 0; r < matResult.rows; r++)
	{
		float* r_matResult = matResult.ptr<float>(r);
		float* r_img_mask = img_mask.ptr<float>(r);
		float* r_Img = Img.ptr<float>(r);
		float* r_templx_mask, * r_sub_source, * r_mask;
		for (int c = 0; c < matResult.cols; ++c, ++r_matResult, ++r_Img, ++r_img_mask)
		{
			r_templx_mask = templx_mask.ptr<float>();
			r_mask = Mask.ptr<float>();
			r_sub_source = r_Img;
			for (t_r = 0; t_r < t_r_end; ++t_r, r_sub_source += Img.cols, r_templx_mask += Templ.cols, r_mask += Templ.cols)
			{
				*r_matResult = *r_matResult + IM_Conv32f_CVSIMD(r_templx_mask, r_sub_source, Templ.cols);
				*r_img_mask = *r_img_mask + IM_Conv32f_CVSIMD(r_mask, r_sub_source, Templ.cols);
			}
		}
	}

	Mat temp_res = img_mask.mul(sum(templx_mask).div(mask_sum));
	matResult -= temp_res;
	double norm_templx = norm(Mask.mul(Templ - sum(Mask.mul(Templ)).div(mask_sum)), NORM_L2);
	Mat norm_imgx(corrSize, CV_32F); norm_imgx.setTo(0);
	Mat img2 = Img.mul(Img);
	Mat mask2 = Mask.mul(Mask);
	Scalar mask2_sum = sum(mask2);
	Mat img_mask2_corr(corrSize, Img.type()); img_mask2_corr.setTo(0);
	for (int r = 0; r < matResult.rows; r++)
	{
		float* r_norm_imgx = norm_imgx.ptr<float>(r);
		float* r_img_mask2_corr = img_mask2_corr.ptr<float>(r);
		float* r_Img = Img.ptr<float>(r);
		float* r_img2 = img2.ptr<float>(r);
		float* r_sub_source1, * r_sub_source2, * r_mask2;
		for (int c = 0; c < matResult.cols; ++c, ++r_norm_imgx, ++r_Img, ++r_img2, ++r_img_mask2_corr)
		{
			r_mask2 = mask2.ptr<float>();
			r_sub_source1 = r_Img;
			r_sub_source2 = r_img2;
			for (t_r = 0; t_r < t_r_end; ++t_r, r_sub_source1 += Img.cols, r_sub_source2 += Img.cols, r_mask2 += Templ.cols)
			{
				*r_norm_imgx = *r_norm_imgx + IM_Conv32f_CVSIMD(r_mask2, r_sub_source2, Templ.cols);
				*r_img_mask2_corr = *r_img_mask2_corr + IM_Conv32f_CVSIMD(r_mask2, r_sub_source1, Templ.cols);
			}
		}
	}

	temp_res = img_mask.mul(Scalar(1.0, 1.0, 1.0, 1.0).div(mask_sum)).mul(img_mask.mul(mask2_sum.div(mask_sum)) - 2 * img_mask2_corr);
	norm_imgx += temp_res;
	cv::sqrt(norm_imgx, norm_imgx);
	matResult = matResult / (norm_imgx * norm_templx);

}

void GrayscaleMatch::MatchTemplate(cv::Mat& matSrc, s_TemplData* pTemplData, cv::Mat& matResult, int iLayer)
{
	matResult.create(matSrc.rows - pTemplData->vecPyramid[iLayer].rows + 1, matSrc.cols - pTemplData->vecPyramid[iLayer].cols + 1, CV_32FC1);
	matResult.setTo(0);
	if (pTemplData->vecResultEqual1[iLayer])
	{
		matResult = Scalar::all(1);
		return;
	}
	cv::Mat& matTemplate = pTemplData->vecPyramid[iLayer];
	Mat sum, sqsum;
	cv::integral(matSrc, sum, sqsum, CV_64F);
	double* p0 = (double*)sum.data;
	double* p1 = p0 + pTemplData->vecPyramid[iLayer].cols;
	double* p2 = (double*)(sum.data + pTemplData->vecPyramid[iLayer].rows * sum.step);
	double* p3 = p2 + pTemplData->vecPyramid[iLayer].cols;
	double* q0 = (double*)sqsum.data;
	double* q1 = q0 + pTemplData->vecPyramid[iLayer].cols;
	double* q2 = (double*)(sqsum.data + pTemplData->vecPyramid[iLayer].rows * sqsum.step);
	double* q3 = q2 + pTemplData->vecPyramid[iLayer].cols;
	int sumstep = sum.data ? (int)(sum.step / sizeof(double)) : 0;
	int sqstep = sqsum.data ? (int)(sqsum.step / sizeof(double)) : 0;
	double dTemplMean0 = pTemplData->vecTemplMean[iLayer][0];
	double dTemplNorm = pTemplData->vecTemplNorm[iLayer];
	double dInvArea = pTemplData->vecInvArea[iLayer];
	for (int r = 0; r < matResult.rows; r++)
	{
		float* r_matResult = matResult.ptr<float>(r);
		uchar* r_source = matSrc.ptr<uchar>(r);
		uchar* r_template, * r_sub_source;
		int idx = r * sumstep;
		int idx2 = r * sqstep;
		if (m_params.useMean || m_params.useSDV)
		{
			for (int c = 0; c < matResult.cols; ++c, ++r_matResult, ++r_source, idx += 1, idx2 += 1)
			{
				bool not_mach = 0, not_mach_mean = 0, not_mach_standard_deviation = 0;
				float t_mean = (p0[idx] - p1[idx] - p2[idx] + p3[idx]) * dInvArea;
				if (m_params.useMean)
				{
					not_mach_mean = (abs(t_mean - dTemplMean0) >= m_params.mean) ? 1 : 0;
				}
				if (m_params.useSDV)
				{
					float standard_deviation = cv::sqrt(abs((q0[idx] - q1[idx] - q2[idx] + q3[idx]) * dInvArea - t_mean * t_mean));
					not_mach_standard_deviation = (abs(standard_deviation - cv::sqrt(dTemplNorm)) >= m_params.SDV) ? 1 : 0;
				}
				not_mach = bool(not_mach_mean + not_mach_standard_deviation);
				if (not_mach_mean || not_mach_standard_deviation)
				{
					continue;
				}
				else
				{
					r_template = matTemplate.ptr<uchar>();
					r_sub_source = r_source;
					for (int t_r = 0; t_r < matTemplate.rows; ++t_r, r_sub_source += matSrc.cols, r_template += matTemplate.cols)
					{
						*r_matResult = *r_matResult + IM_Conv_CVSIMD(r_template, r_sub_source, matTemplate.cols);
					}
					double num = *r_matResult, t;
					double wndMean2 = 0, wndSum2 = 0;
					t = p0[idx] - p1[idx] - p2[idx] + p3[idx];
					wndMean2 += t * t;
					num -= t * dTemplMean0;
					wndMean2 *= dInvArea;

					t = q0[idx2] - q1[idx2] - q2[idx2] + q3[idx2];
					wndSum2 += t;
					double diff2 = MAX(wndSum2 - wndMean2, 0);
					if (diff2 <= (std::min)(0.5, 10 * FLT_EPSILON * wndSum2)) { t = 0; } // avoid rounding errors
					else { t = std::sqrt(diff2) * dTemplNorm; }
					if (fabs(num) < t) { num /= t; }
					else if (fabs(num) < t * 1.125) { num = num > 0 ? 1 : -1; }
					else { num = 0; }
					*r_matResult = (float)num;
				}
			}
		}
		else
		{
			for (int c = 0; c < matResult.cols; ++c, ++r_matResult, ++r_source, idx += 1, idx2 += 1)
			{
				r_template = matTemplate.ptr<uchar>();
				r_sub_source = r_source;
				for (int t_r = 0; t_r < matTemplate.rows; ++t_r, r_sub_source += matSrc.cols, r_template += matTemplate.cols)
				{
					*r_matResult = *r_matResult + IM_Conv_CVSIMD(r_template, r_sub_source, matTemplate.cols);
				}
				double num = *r_matResult, t;
				double wndMean2 = 0, wndSum2 = 0;
				t = p0[idx] - p1[idx] - p2[idx] + p3[idx];
				wndMean2 += t * t;
				num -= t * dTemplMean0;
				wndMean2 *= dInvArea;

				t = q0[idx2] - q1[idx2] - q2[idx2] + q3[idx2];
				wndSum2 += t;
				double diff2 = MAX(wndSum2 - wndMean2, 0);
				if (diff2 <= (std::min)(0.5, 10 * FLT_EPSILON * wndSum2)) { t = 0; } // avoid rounding errors
				else { t = std::sqrt(diff2) * dTemplNorm; }
				if (fabs(num) < t) { num /= t; }
				else if (fabs(num) < t * 1.125) { num = num > 0 ? 1 : -1; }
				else { num = 0; }
				*r_matResult = (float)num;
			}
		}
	}

}

// 顶层匹配


void GrayscaleMatch::topMatchThread(std::vector<float> vecAngles, int start, int end, cv::Point2f ptCenter,
	std::vector<cv::Mat1b>& vecMatSrcPyr, s_TemplData& pTemplData,
	std::vector<float> vecLayerScore, std::vector<s_MatchParameter>& vecMatchParameter)
{
	for (int i = start; i < end; i++)
	{
		Mat matRotatedSrc, matR = getRotationMatrix2D(ptCenter, vecAngles[i], 1);
		Mat1f matResult;
		Point ptMaxLoc;
		double dValue, dMaxVal;
		Size sizeBest = GetBestRotationSize(vecMatSrcPyr[m_params.pyramidLayer].size(), pTemplData.vecPyramid[m_params.pyramidLayer].size(), vecAngles[i]);

		float fTranslationX = (sizeBest.width - 1) / 2.0f - ptCenter.x;
		float fTranslationY = (sizeBest.height - 1) / 2.0f - ptCenter.y;
		matR.at<double>(0, 2) += fTranslationX;
		matR.at<double>(1, 2) += fTranslationY;
		warpAffine(vecMatSrcPyr[m_params.pyramidLayer], matRotatedSrc, matR, sizeBest, INTER_LINEAR, BORDER_CONSTANT, Scalar(pTemplData.iBorderColor));
		if (m_params.useMask)
		{
			MatchTemplateMask(matRotatedSrc, &pTemplData, matResult, m_params.pyramidLayer);
		}
		else
		{
			MatchTemplate(matRotatedSrc, &pTemplData, matResult, m_params.pyramidLayer);
		}

		int size = 3;
		for (int r = 0; r < matResult.rows; r++)
		{
			for (int c = 0; c < matResult.cols; c++)
			{
				float score = matResult[r][c];
				bool flag = false;
				if (score > vecLayerScore[m_params.pyramidLayer])
				{
					float maxScore = 0;
					for (int m = max(0, r - size); m < min(matResult.rows, r + size + 1); m++)
					{
						for (int n = max(0, c - size); n < min(matResult.cols, c + size + 1); n++)
						{
							maxScore = maxScore > matResult[m][n] ? maxScore : matResult[m][n];
						}
					}
					if (score >= maxScore)
					{
						vecMatchParameter.push_back(s_MatchParameter(Point2f(c - fTranslationX, r - fTranslationY), score, vecAngles[i]));
					}
				}
			}
		}
	}

	

}

void GrayscaleMatch::topMatch(vector<s_MatchParameter>& _topMatchParameters)
{
	int iTopLayer = m_params.pyramidLayer;
	// 匹配角度
	std::vector<float> t_topAngs;
	if (m_params.angleRange > 0)
	{
		for (double dAngle = 0; dAngle < m_params.angleRange; dAngle += m_anglePyr[iTopLayer])
			t_topAngs.push_back(m_params.angleStart + dAngle);
	}
	else
	{
		t_topAngs.push_back(m_params.angleStart);
	}

	int iTopSrcW = m_grayPyr[iTopLayer].cols, iTopSrcH = m_grayPyr[iTopLayer].rows;
	Point2f ptCenter((iTopSrcW - 1) / 2.0f, (iTopSrcH - 1) / 2.0f);
	int iSize = (int)t_topAngs.size();

	unsigned int const hardware_threads = std::thread::hardware_concurrency();
	unsigned int const min_per_thread = 1;
	unsigned int const max_threads = (iSize + min_per_thread - 1) / min_per_thread;
	unsigned int const num_threads = std::min(hardware_threads != 0 ? hardware_threads : 4, max_threads);
	unsigned int const block_size = iSize / num_threads;
	std::vector<vector<s_MatchParameter>> vecMatchParameter_s(num_threads);
	std::vector<std::thread> threads(num_threads - 1);
	int start = 0;
	int end = 0;
	for (unsigned long i = 0; i < (num_threads - 1); ++i)
	{
		end = start + block_size;
		threads[i] = thread(&GrayscaleMatch::topMatchThread, this, t_topAngs, start, end, 
			ptCenter, std::ref(m_grayPyr), std::ref(m_TemplData), m_scorePyr, std::ref(vecMatchParameter_s[i]));
		start = end;
	}
	topMatchThread(t_topAngs, start, iSize, ptCenter, m_grayPyr, m_TemplData, m_scorePyr, vecMatchParameter_s[num_threads - 1]);
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

	_topMatchParameters.clear();
	for (int i = 0; i < num_threads; i++)
	{
		vector<s_MatchParameter> temp = vecMatchParameter_s[i];
		_topMatchParameters.insert(_topMatchParameters.end(), temp.begin(), temp.end());
	}

	int iMatchSize = (int)_topMatchParameters.size();
	int iDstW = m_TemplData.vecPyramid[m_params.pyramidLayer].cols, iDstH = m_TemplData.vecPyramid[m_params.pyramidLayer].rows;
	for (int i = 0; i < iMatchSize; i++)
	{
		Point2f ptLT, ptRT, ptRB, ptLB;
		double dRAngle = -_topMatchParameters[i].dMatchAngle * D2R;
		ptLT = ptRotatePt2f(_topMatchParameters[i].pt, ptCenter, dRAngle);
		ptRT = Point2f(ptLT.x + iDstW * (float)cos(dRAngle), ptLT.y - iDstW * (float)sin(dRAngle));
		ptLB = Point2f(ptLT.x + iDstH * (float)sin(dRAngle), ptLT.y + iDstH * (float)cos(dRAngle));
		ptRB = Point2f(ptRT.x + iDstH * (float)sin(dRAngle), ptRT.y + iDstH * (float)cos(dRAngle));
		Point2f ptRectCenter = Point2f((ptLT.x + ptRT.x + ptLB.x + ptRB.x) / 4.0f, (ptLT.y + ptRT.y + ptLB.y + ptRB.y) / 4.0f);
		_topMatchParameters[i].rectR = RotatedRect(ptRectCenter, m_TemplData.vecPyramid[m_params.pyramidLayer].size(), (float)_topMatchParameters[i].dMatchAngle);
	}
}


void GetRotatedROI(Mat& matSrc, Size size, Point2f ptLT, double dAngle, Mat& matROI)
{
	double dAngle_radian = dAngle * D2R;
	Point2f ptC((matSrc.cols - 1) / 2.0f, (matSrc.rows - 1) / 2.0f);
	Point2f ptLT_rotate = ptRotatePt2f(ptLT, ptC, dAngle_radian);
	Size sizePadding(size.width + 6, size.height + 6);


	Mat rMat = getRotationMatrix2D(ptC, dAngle, 1);
	rMat.at<double>(0, 2) -= ptLT_rotate.x - 3;
	rMat.at<double>(1, 2) -= ptLT_rotate.y - 3;

	warpAffine(matSrc, matROI, rMat, sizePadding);
}

void GrayscaleMatch::bottomMatchThread(std::vector<s_MatchParameter>& vecMatchParameter, int start, int end, Point2f ptCenter,
	int iDstW, int iDstH, int iStopLayer, std::vector<cv::Mat1b>& vecMatSrcPyr, s_TemplData& pTemplData,
	 std::vector<float> vecLayerScore, vector<s_MatchParameter>& vecAllResult)
{

	for (int i = start; i < end; i++)
	{
		double dRAngle = -vecMatchParameter[i].dMatchAngle * D2R;
		Point2f ptLT = ptRotatePt2f(vecMatchParameter[i].pt, ptCenter, dRAngle);
		double dAngleStep = atan(1.0 / max(iDstW, iDstH)) * R2D;
		vecMatchParameter[i].dAngleStart = vecMatchParameter[i].dMatchAngle - dAngleStep;
		vecMatchParameter[i].dAngleEnd = vecMatchParameter[i].dMatchAngle + dAngleStep;
		if (m_params.pyramidLayer <= iStopLayer)
		{
			vecMatchParameter[i].pt = Point2d(ptLT * ((m_params.pyramidLayer == 0) ? 1 : 2));
			vecAllResult.push_back(vecMatchParameter[i]);
		}
		else
		{
			for (int iLayer = m_params.pyramidLayer - 1; iLayer >= iStopLayer; iLayer--)
			{
				dAngleStep = atan(1.0 / max(pTemplData.vecPyramid[iLayer].cols, pTemplData.vecPyramid[iLayer].rows)) * R2D;
				vector<double> vecAngles;
				double dMatchedAngle = vecMatchParameter[i].dMatchAngle;
				for (int i = -1; i <= 1; i++)
					vecAngles.push_back(dMatchedAngle + dAngleStep * i);

				Point2f ptSrcCenter((vecMatSrcPyr[iLayer].cols - 1) / 2.0f, (vecMatSrcPyr[iLayer].rows - 1) / 2.0f);
				int iSize = (int)vecAngles.size();
				vector<s_MatchParameter> vecNewMatchParameter(iSize);
				int iMaxScoreIndex = 0;
				double dBigValue = -1;
				for (int j = 0; j < iSize; j++)
				{
					Mat matResult, matRotatedSrc;
					double dMaxValue = 0;
					Point ptMaxLoc;
					GetRotatedROI(vecMatSrcPyr[iLayer], pTemplData.vecPyramid[iLayer].size(), ptLT * 2, vecAngles[j], matRotatedSrc);
					if (m_params.useMask)
					{
						MatchTemplateMask(matRotatedSrc, &pTemplData, matResult, iLayer);
					}
					else
					{
						MatchTemplate(matRotatedSrc, &pTemplData, matResult, iLayer);
					}
					minMaxLoc(matResult, 0, &dMaxValue, 0, &ptMaxLoc);
					vecNewMatchParameter[j] = s_MatchParameter(ptMaxLoc, dMaxValue, vecAngles[j]);

					if (vecNewMatchParameter[j].dMatchScore > dBigValue)
					{
						iMaxScoreIndex = j;
						dBigValue = vecNewMatchParameter[j].dMatchScore;
					}
					if (ptMaxLoc.x == 0 || ptMaxLoc.y == 0 || ptMaxLoc.y == matResult.cols - 1 || ptMaxLoc.x == matResult.rows - 1)
						vecNewMatchParameter[j].bPosOnBorder = true;
					if (!vecNewMatchParameter[j].bPosOnBorder)
					{
						for (int y = -1; y <= 1; y++)
							for (int x = -1; x <= 1; x++)
								vecNewMatchParameter[j].vecResult[x + 1][y + 1] = matResult.at<float>(ptMaxLoc + Point(x, y));
					}
				}

				if (vecNewMatchParameter[iMaxScoreIndex].dMatchScore < vecLayerScore[iLayer])
					break;
	
				double dNewMatchAngle = vecNewMatchParameter[iMaxScoreIndex].dMatchAngle;
				Point2f ptPaddingLT = ptRotatePt2f(ptLT * 2, ptSrcCenter, dNewMatchAngle * D2R) - Point2f(3, 3);
				Point2f pt(vecNewMatchParameter[iMaxScoreIndex].pt.x + ptPaddingLT.x, vecNewMatchParameter[iMaxScoreIndex].pt.y + ptPaddingLT.y);
				pt = ptRotatePt2f(pt, ptSrcCenter, -dNewMatchAngle * D2R);
				if (iLayer == iStopLayer)
				{
					vecNewMatchParameter[iMaxScoreIndex].pt = pt * (iStopLayer == 0 ? 1 : 2);
					vecAllResult.push_back(vecNewMatchParameter[iMaxScoreIndex]);
				}
				else
				{
					vecMatchParameter[i].dMatchAngle = dNewMatchAngle;
					vecMatchParameter[i].dAngleStart = vecMatchParameter[i].dMatchAngle - dAngleStep / 2;
					vecMatchParameter[i].dAngleEnd = vecMatchParameter[i].dMatchAngle + dAngleStep / 2;
					ptLT = pt;
				}
			}

		}
	}
}

void GrayscaleMatch::bottomMatch(std::vector<s_MatchParameter>& in_MatchParameters, std::vector<s_MatchParameter>& out_MatchParameters)
{
	out_MatchParameters.clear();
	int iMatchSize = in_MatchParameters.size();
	int iTopSrcW = m_grayPyr[m_params.pyramidLayer].cols, iTopSrcH = m_grayPyr[m_params.pyramidLayer].rows;
	int iDstW = m_TemplData.vecPyramid[m_params.pyramidLayer].cols, iDstH = m_TemplData.vecPyramid[m_params.pyramidLayer].rows;

	Point2f ptCenter((iTopSrcW - 1) / 2.0f, (iTopSrcH - 1) / 2.0f);
	unsigned int const hardware_threads = std::thread::hardware_concurrency();
	unsigned int const min_per_thread = 1;
	unsigned int const max_threads = (iMatchSize + min_per_thread - 1) / min_per_thread;
	unsigned int const num_threads = std::min(hardware_threads != 0 ? hardware_threads : 4, max_threads);
	unsigned int const block_size = iMatchSize / num_threads;
	std::vector<vector<s_MatchParameter>> vecAllResult_s(num_threads);
	std::vector<std::thread> threads(num_threads - 1);
	int start = 0;
	int end = 0;
	int iStopLayer = 0;
	for (unsigned long i = 0; i < (num_threads - 1); ++i)
	{
		end = start + block_size;
		threads[i] = thread(&GrayscaleMatch::bottomMatchThread, this, std::ref(in_MatchParameters), start, end, ptCenter
			, iDstW, iDstH, iStopLayer, std::ref(m_grayPyr), std::ref(m_TemplData), m_scorePyr
			, std::ref(vecAllResult_s[i]));
		start = end;
	}

	bottomMatchThread(std::ref(in_MatchParameters), start, iMatchSize, ptCenter, iDstW, iDstH, iStopLayer, m_grayPyr,
		m_TemplData, m_scorePyr, std::ref(vecAllResult_s[num_threads - 1]));
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
	for (int i = 0; i < num_threads; i++)
	{
		vector<s_MatchParameter> temp = vecAllResult_s[i];
		out_MatchParameters.insert(out_MatchParameters.end(), temp.begin(), temp.end());
	}
}


void FilterWithScore(vector<s_MatchParameter>* vec, float dScore)
{
	std::sort(vec->begin(), vec->end(), [](const s_MatchParameter& lhs, const s_MatchParameter& rhs) { return  lhs.dMatchScore > rhs.dMatchScore; });
	int iSize = vec->size(), iIndexDelete = iSize + 1;
	for (int i = 0; i < iSize; i++)
	{
		if ((*vec)[i].dMatchScore < dScore)
		{
			iIndexDelete = i;
			break;
		}
	}
	if (iIndexDelete == iSize + 1)
		return;
	vec->erase(vec->begin() + iIndexDelete, vec->end());
	return;
	vector<s_MatchParameter>::iterator it;
	for (it = vec->begin(); it != vec->end();)
	{
		if (((*it).dMatchScore < dScore))
			it = vec->erase(it);
		else
			++it;
	}
}

void FilterWithRotatedRect(vector<s_MatchParameter>& vec, int iMethod, double dMaxOverLap)
{
	int iMatchSize = (int)vec.size();
	RotatedRect rect1, rect2;
	for (int i = 0; i < iMatchSize - 1; i++)
	{
		if (vec[i].bDelete)
			continue;

		for (int j = i + 1; j < iMatchSize; j++)
		{
			if (vec[j].bDelete)
				continue;
			rect1 = vec[i].rectR;
			rect2 = vec[j].rectR;
			vector<Point2f> vecInterSec;
			int iInterSecType = rotatedRectangleIntersection(rect1, rect2, vecInterSec);
			if (iInterSecType == INTERSECT_NONE)
				continue;
			else if (iInterSecType == INTERSECT_FULL)
			{
				int iDeleteIndex;
				/*if (iMethod == TM_SQDIFF)
					iDeleteIndex = (vec[i].dMatchScore <= vec[j].dMatchScore) ? j : i;
				else*/
				iDeleteIndex = (vec[i].dMatchScore >= vec[j].dMatchScore) ? j : i;
				vec[iDeleteIndex].bDelete = true;
			}
			else
			{
				if (vecInterSec.size() < 3)
					continue;
				else
				{
					int iDeleteIndex;
					//SortPtWithCenter(vecInterSec);
					double dArea = contourArea(vecInterSec);
					double dRatio = dArea / rect1.size.area();
					if (dRatio > dMaxOverLap)
					{
						/*if (iMethod == TM_SQDIFF)
							iDeleteIndex = (vec[i].dMatchScore <= vec[j].dMatchScore) ? j : i;
						else*/
						iDeleteIndex = (vec[i].dMatchScore >= vec[j].dMatchScore) ? j : i;
						vec[iDeleteIndex].bDelete = true;
					}
				}
			}
		}
	}
	vector<s_MatchParameter>::iterator it;
	for (it = vec.begin(); it != vec.end();)
	{
		if ((*it).bDelete)
			it = vec.erase(it);
		else
			++it;
	}
}

void GrayscaleMatch::filtResult(std::vector<s_MatchParameter>& in_MatchParameters, std::vector<s_SingleTargetMatch>& out_MatchParameters)
{
	FilterWithScore(&in_MatchParameters, m_params.score);

	int iDstW = m_TemplData.vecPyramid[0].cols, iDstH = m_TemplData.vecPyramid[0].rows;
	for (int i = 0; i < (int)in_MatchParameters.size(); i++)
	{
		Point2f ptLT, ptRT, ptRB, ptLB;
		double dRAngle = -in_MatchParameters[i].dMatchAngle * D2R;
		ptLT = in_MatchParameters[i].pt;
		ptRT = Point2f(ptLT.x + iDstW * (float)cos(dRAngle), ptLT.y - iDstW * (float)sin(dRAngle));
		ptLB = Point2f(ptLT.x + iDstH * (float)sin(dRAngle), ptLT.y + iDstH * (float)cos(dRAngle));
		ptRB = Point2f(ptRT.x + iDstH * (float)sin(dRAngle), ptRT.y + iDstH * (float)cos(dRAngle));
		Point2f ptRectCenter = Point2f((ptLT.x + ptRT.x + ptLB.x + ptRB.x) / 4.0f, (ptLT.y + ptRT.y + ptLB.y + ptRB.y) / 4.0f);
		in_MatchParameters[i].rectR = RotatedRect(ptRectCenter, m_TemplData.vecPyramid[0].size(), (float)in_MatchParameters[i].dMatchAngle);
	}
	double FilterWithRotatedRect_start = getTickCount();
	FilterWithRotatedRect(in_MatchParameters, TM_CCOEFF_NORMED, m_params.maxOverlap);
	std::sort(in_MatchParameters.begin(), in_MatchParameters.end(), [](const s_MatchParameter& lhs, const s_MatchParameter& rhs) { return  lhs.dMatchScore > rhs.dMatchScore; });
	out_MatchParameters.clear();
	int iMatchSize = (int)in_MatchParameters.size();
	if (in_MatchParameters.size() == 0)
		return;
	int iW = m_TemplData.vecPyramid[0].cols, iH = m_TemplData.vecPyramid[0].rows;
	for (int i = 0; i < iMatchSize; i++)
	{
		s_SingleTargetMatch sstm;
		double dRAngle = -in_MatchParameters[i].dMatchAngle * D2R;
		sstm.ptLT = in_MatchParameters[i].pt;
		sstm.ptRT = Point2d(sstm.ptLT.x + iW * cos(dRAngle), sstm.ptLT.y - iW * sin(dRAngle));
		sstm.ptLB = Point2d(sstm.ptLT.x + iH * sin(dRAngle), sstm.ptLT.y + iH * cos(dRAngle));
		sstm.ptRB = Point2d(sstm.ptRT.x + iH * sin(dRAngle), sstm.ptRT.y + iH * cos(dRAngle));
		sstm.ptCenter = Point2d((sstm.ptLT.x + sstm.ptRT.x + sstm.ptRB.x + sstm.ptLB.x) / 4, (sstm.ptLT.y + sstm.ptRT.y + sstm.ptRB.y + sstm.ptLB.y) / 4);
		sstm.dMatchedAngle = -in_MatchParameters[i].dMatchAngle;
		sstm.dMatchScore = in_MatchParameters[i].dMatchScore * 100;
		/*0 -- 360*/
		if (sstm.dMatchedAngle < -180)
			sstm.dMatchedAngle += 360;
		if (sstm.dMatchedAngle > 180)
			sstm.dMatchedAngle -= 360;
		/*if (sstm.dMatchedAngle < 0)
			sstm.dMatchedAngle += 360;*/
		out_MatchParameters.push_back(sstm);
		if (i + 1 == m_params.maxtargs)
			break;
	}
}


// 匹配
void GrayscaleMatch::match(cv::Mat _image)
{
	// 0 参数初始化
	m_result = Mat(1, 1, CV_32FC1, Scalar(0));
	m_gray = _image.clone();
	bool t_status = initialize();
	if (!t_status)
	{
		return;
	}

	// 1 顶层匹配
	vector<s_MatchParameter> t_topMatchParameters;
	topMatch(t_topMatchParameters);
	if (!t_topMatchParameters.size())
	{
		return;
	}

	// 2 其他层匹配
	vector<s_MatchParameter> t_bottomMatchParameters;
	bottomMatch(t_topMatchParameters, t_bottomMatchParameters);
	if (!t_bottomMatchParameters.size())
	{
		return;
	}

	// 3 结果处理
	std::vector<s_SingleTargetMatch> t_vecSingleTargetData;;
	filtResult(t_bottomMatchParameters, t_vecSingleTargetData);
	if (!t_vecSingleTargetData.size())
	{
		return;
	}
	Mat1f result1 = Mat(t_vecSingleTargetData.size(), 12, CV_32FC1, Scalar(0));
	for (int i = 0; i < t_vecSingleTargetData.size(); i++)
	{
		cout << t_vecSingleTargetData[i].ptCenter << endl;
		result1[i][0] = ((float)(t_vecSingleTargetData[i].ptCenter.x));
		result1[i][1] = ((float)(t_vecSingleTargetData[i].ptCenter.y));

		result1[i][2] = ((float)(t_vecSingleTargetData[i].dMatchedAngle));
		result1[i][3] = ((float)(t_vecSingleTargetData[i].dMatchScore));


		result1[i][4] = (float)(t_vecSingleTargetData[i].ptLT.x);
		result1[i][5] = ((float)(t_vecSingleTargetData[i].ptLT.y));

		result1[i][6] = ((float)(t_vecSingleTargetData[i].ptLB.x));
		result1[i][7] = ((float)(t_vecSingleTargetData[i].ptLB.y));

		result1[i][8] = ((float)(t_vecSingleTargetData[i].ptRB.x));
		result1[i][9] = ((float)(t_vecSingleTargetData[i].ptRB.y));

		result1[i][10] = ((float)(t_vecSingleTargetData[i].ptRT.x));
		result1[i][11] = ((float)(t_vecSingleTargetData[i].ptRT.y));
	}
	m_result = result1;

#if 1
	cout << t_topMatchParameters.size() << endl;
	cout << t_bottomMatchParameters.size() << endl;
	cout << t_vecSingleTargetData.size() << endl;
#endif

}