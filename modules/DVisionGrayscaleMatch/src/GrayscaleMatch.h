#pragma once

#include <opencv2/opencv.hpp>

struct Params
{
	Params() {};
	cv::Mat temp;				// 模板图像

	bool useMask = 0;			// 使用掩码
	cv::Mat mask = cv::Mat();	// 掩码

	float score = 0.8;			// 匹配分数

	int maxtargs = 0;			// 最大匹配数量
	int pyramidLayer = 256;		// 金字塔参数

	float angleStart = 0;		// 匹配角度起始值
	float angleRange = 360;		// 匹配角度范围
	float maxOverlap = 0;		// 最大重叠率

	bool useMean = 0;			// 均值加速
	float mean = 30;
	bool useSDV = 0;			// 标准差加速
	float SDV = 30;

	bool polarity = true;		// 极性
	bool subpixel = 0;			// 亚像素
};

struct s_TemplData
{
	std::vector<cv::Mat> vecPyramid;				/*Template image Pyramid*/
	std::vector<cv::Scalar> vecTemplMean;		/*Template pyramid image average*/
	std::vector<double> vecTemplNorm;		/*Template pyramid image standard deviation*/
	std::vector<double> vecInvArea;			/* 1 / (Template pyramid image area) */
	std::vector<bool> vecResultEqual1;
	bool bIsPatternLearned;
	int iBorderColor;
	void clear()
	{
		std::vector<cv::Mat>().swap(vecPyramid);
		std::vector<double>().swap(vecTemplNorm);
		std::vector<double>().swap(vecInvArea);
		std::vector<cv::Scalar>().swap(vecTemplMean);
		std::vector<bool>().swap(vecResultEqual1);
	}
	void resize(int iSize)
	{
		vecTemplMean.resize(iSize);
		vecTemplNorm.resize(iSize, 0);
		vecInvArea.resize(iSize, 1);
		vecResultEqual1.resize(iSize, false);
	}
	s_TemplData()
	{
		bIsPatternLearned = false;
	}
};

/*Match Parameter struct*/
struct s_MatchParameter
{
	cv::Point2d pt;					/*Match result point*/
	double dMatchScore;			/*Match result Score*/
	double dMatchAngle;			/*Match result Angle*/
	double dAngleStart;			/*Match result Angle range*/
	double dAngleEnd;
	cv::RotatedRect rectR;			/*OpenCV RotatedRect*/
	bool bDelete;				/*do match?*/

	double vecResult[3][3];		/*for subpixel*/
	int iMaxScoreIndex;
	bool bPosOnBorder;
	cv::Point2d ptSubPixel;
	double dNewAngle;

	s_MatchParameter(cv::Point2f ptMinMax, double dScore, double dAngle)//, Mat matRotatedSrc = Mat ())
	{
		pt = ptMinMax;
		dMatchScore = dScore;
		dMatchAngle = dAngle;

		bDelete = false;
		dNewAngle = 0.0;

		bPosOnBorder = false;
	}
	s_MatchParameter()
	{
		double dMatchScore = 0;
		double dMatchAngle = 0;
	}
	~s_MatchParameter()
	{

	}
};

struct s_SingleTargetMatch
{
	cv::Point2d ptLT, ptRT, ptRB, ptLB, ptCenter;
	double dMatchedAngle;
	double dMatchScore;
};

class GrayscaleMatch
{
public:
	GrayscaleMatch() {}
	~GrayscaleMatch() {}

public:
	void setTempimage(cv::Mat _temp, bool useMask = 0, cv::Mat _mask = cv::Mat());
	void setScore(float score);
	void setMaxtargs(float maxtargs);
	void setPyramidLayer(float pyramidLayer);
	void setAngle(float angleStart, float angleRange);
	void setMaxOverlap(float maxOverlap);
	void setPolarity(bool polarity);
	void setSubpixel(bool subpixel);
	void setMean(bool useMean, float mean);
	void setSDV(bool useSDV, float SDV);

public:
	void match(cv::Mat _image);
	cv::Mat getResult() { return m_result; }

private:
	bool initialize();

	void MatchTemplateMask(cv::Mat& matSrc, s_TemplData* pTemplData, cv::Mat& matResult, int iLayer);
	void MatchTemplate(cv::Mat& matSrc, s_TemplData* pTemplData, cv::Mat& matResult, int iLayer);

	void topMatchThread(std::vector<float> vecAngles, int start, int end, cv::Point2f ptCenter,
						std::vector<cv::Mat1b>& vecMatSrcPyr, s_TemplData& pTemplData, 
						std::vector<float> vecLayerScore, std::vector<s_MatchParameter>& vecMatchParameter);
	void topMatch(std::vector<s_MatchParameter>& _topMatchParameters);

	void bottomMatchThread(std::vector<s_MatchParameter>& vecMatchParameter, int start, int end, cv::Point2f ptCenter, 
		int iDstW, int iDstH, int iStopLayer, std::vector<cv::Mat1b>& vecMatSrcPyr, s_TemplData& pTemplData, 
		std::vector<float> vecLayerScore, std::vector<s_MatchParameter>& vecAllResult);
	void bottomMatch(std::vector<s_MatchParameter>& in_MatchParameters, std::vector<s_MatchParameter>& out_MatchParameters);

	void filtResult(std::vector<s_MatchParameter>& in_MatchParameters, std::vector<s_SingleTargetMatch>& out_MatchParameters);

private:
	float m_angleStep = 1;
	cv::Mat m_gray;
	cv::Mat1f m_result;
	std::vector<float> m_anglePyr;
	std::vector<float> m_scorePyr;
	std::vector<cv::Mat1b> m_maskPyr;
	std::vector<cv::Mat1b> m_grayPyr;

	Params m_params = Params();
	s_TemplData m_TemplData;



};