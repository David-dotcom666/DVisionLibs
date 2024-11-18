#pragma once

#include<opencv2/opencv.hpp>


struct Params
{
	float inlier = 0.5;					// 弧段的内点率
	float closedInlier = 0.5;			// 闭合圆的内点率
	float lineLen = 160;					// 直线弧线交界处断开
	float centerDis = 5;				// 圆心距离
	float radiusDis = 5;				// 半径差距
	float sharpAngle = 60;				// 角度急转处断开
	float SegmentCurvature = 0.01;	// 弧线段曲率，根据弯曲程度断开
	float minLen = 32;					// 最小边缘段
	int filtSize = 5;					// 滤波尺寸
};

struct CircleData
{
	float radius = 0;							// 半径
	cv::Point2f center;							// 圆心
	float ratio = 1;							// 内点率

};

class CircleDetect
{
public:
	CircleDetect() {};
	~CircleDetect() {};

public:
	void setSplit(float lineLen, float sharpAngle, float SegmentCurvature);
	void setInlier(float inlier, float closedInlier);
	void setDistance(float centerDis, float radiusDis);
	void setFiltSize(int filtSize);
	void setMinLen(float minLen);

public:
	void detect(cv::Mat image);
	cv::Mat getResult() {return m_results;}

private:
	void extractClosedEdges(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& out_sigments);
	void RDP(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& out_sigments);
	void rejectSharpTurn(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
						 std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists);
	void detectline(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
					std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists);
	void detectInflexPt(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
						std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists);
	void detectshort(std::vector<std::vector<cv::Point>>& in_out_edgelists);
	void extractClosedEdges2(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& out_edgelists
							, std::vector<std::vector<cv::Point>>& out_closededgelists);
	void coCircleGroupArcs(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& out_groupedArcs,
		std::vector<std::vector<cv::Point>>& out_groupedArcsThreePt, std::vector<cv::Vec3f>& out_recordOR);
	void toalcircle(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& in_groupedArcs, std::vector<std::vector<cv::Point>>& in_groupedArcsThreePt,
		std::vector<cv::Vec3f>& in_recordOR, std::vector<CircleData>& out_Circles);
	void clusterCircles(std::vector<CircleData>& in_Circles, std::vector<CircleData>& out_Circles);

private:
	cv::Mat m_gray;
	cv::Mat1f m_results;
	Params m_params = Params();

};
