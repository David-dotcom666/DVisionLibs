#include "BlobDetector.h"
#include <opencv2/core/simd_intrinsics.hpp>

using namespace cv;
using namespace std;


void BlobDetector::setColor(bool blobColor)
{
	m_params.blobColor = blobColor;
}

void BlobDetector::setmaxTargets(int maxTargets)
{
	m_params.maxTargets = maxTargets;
}

void BlobDetector::setSort(int sortMethod, int sortFeature)
{
	m_params.sortMethod = sortMethod;
	m_params.sortFeature = sortFeature;
}

void BlobDetector::setThreshold(int autoThreshold, int minThreshold, int maxThreshold)
{
	m_params.autoThreshold = autoThreshold;
	m_params.maxThreshold = maxThreshold;
	m_params.minThreshold = minThreshold;
}

void BlobDetector::setArea(bool filterByArea, int minArea, int maxArea)
{
	m_params.filterByArea = filterByArea;
	m_params.maxArea = maxArea;
	m_params.minArea = minArea;
}

void BlobDetector::setCircularity(bool filterByCircularity, float minCircularity, float maxCircularity)
{
	m_params.filterByCircularity = filterByCircularity;
	m_params.minCircularity = minCircularity;
	m_params.maxCircularity = maxCircularity;
}

void BlobDetector::setInertiaRatio(bool filterByInertia, float minInertiaRatio, float maxInertiaRatio)
{
	m_params.filterByInertia = filterByInertia;
	m_params.minInertiaRatio = minInertiaRatio;
	m_params.maxInertiaRatio = maxInertiaRatio;
}

void BlobDetector::setRectangularity(bool filterByRectangularity, float minRectangularity, float maxRectangularity)
{
	m_params.filterByRectangularity = filterByRectangularity;
	m_params.minRectangularity = minRectangularity;
	m_params.maxRectangularity = maxRectangularity;
}

void BlobDetector::setPerimeter(bool filterByPerimeter, int minPerimeter, int maxPerimeter)
{
	m_params.filterByRectangularity = filterByPerimeter;
	m_params.minRectangularity = minPerimeter;
	m_params.maxRectangularity = maxPerimeter;
}

void BlobDetector::setCentroidOffset(bool filterByCentroidOffset, int minCentroidOffset, int maxCentroidOffset)
{
	m_params.filterByCentroidOffset = filterByCentroidOffset;
	m_params.minCentroidOffset = minCentroidOffset;
	m_params.maxCentroidOffset = maxCentroidOffset;
}



void BlobDetector::detect(cv::Mat _image, cv::Mat _mask)
{
	//0 图像灰度化
	cv::Mat m_gray;
	if (_image.empty())
	{
		std::cout << "image is empty!" << endl;
		return;
	}
	switch (_image.channels())
	{
	case 1:
		m_gray = _image;
		break;
	case 3:
		cvtColor(_image, m_gray, cv::COLOR_BGR2GRAY);
		break;
	case 4:
		cvtColor(_image, m_gray, cv::COLOR_BGRA2GRAY);
		break;
	default:
		break;
	}

	//1 图像二值化
	cv::Mat m_binary = Mat(m_gray.size(), CV_8UC1, Scalar(0));
	if (m_params.autoThreshold == 0)
	{
		if (m_params.maxThreshold == 255)
		{
			threshold(m_gray, m_binary, m_params.minThreshold, 255, THRESH_BINARY);
		}
		else
		{
			if (m_params.maxThreshold >= m_params.minThreshold)
			{
				parallel_for_(Range(0, m_gray.rows), [&](const Range& r)
					{
						for (int row = r.start; row < r.end; row++)
						{
							uchar* gptr = m_gray.ptr<uchar>(row);
							uchar* bptr = m_binary.ptr<uchar>(row);
							for (int col = 0; col < m_gray.cols; col ++)
							{
								if ((gptr[col] > m_params.minThreshold) && (gptr[col] < m_params.maxThreshold))
									bptr[col] = 255;
							}
						}
					});
			}
			else
			{
				parallel_for_(Range(0, m_gray.rows), [&](const Range& r)
					{
						for (int row = r.start; row < r.end; row++)
						{
							uchar* gptr = m_gray.ptr<uchar>(row);
							uchar* bptr = m_binary.ptr<uchar>(row);
							for (int col = 0; col < m_gray.cols; col++)
							{
								if ((gptr[col] > m_params.minThreshold) || (gptr[col] < m_params.maxThreshold))
									bptr[col] = 255;
							}
						}
					});
			}
		}
	}
	else if (m_params.autoThreshold == 2)
	{
		threshold(m_gray, m_binary, 0, 255, THRESH_TRIANGLE);
	}
	else
	{
		threshold(m_gray, m_binary, 0, 255, THRESH_OTSU);
	}
	
	if (!m_params.blobColor)
	{
		m_binary = 255 - m_binary;
	}

	m_binaryImage = m_binary.clone();

	if (!_mask.empty())
	{
		threshold(_mask, _mask, 0, 1, THRESH_BINARY);
		m_binary = m_binary.mul(_mask);
	}

	// 2 获取连通域
	cv::Mat m_connect;
	int connectNumber = connectedComponents(m_binary, m_connect, 8, CV_32S);
	vector<vector<Point>> points(connectNumber-1);
	for (int i = 0; i < points.size(); i++)
	{
		points[i].reserve(m_binary.rows * m_binary.cols * 0.1);
	}
	for (int row = 0; row < m_connect.rows; ++row)
	{
		uint32_t* ptr = m_connect.ptr<uint32_t>(row);
		for (int col = 0; col < m_connect.cols; ++col)
		{
			int label = ptr[col];
			if (label > 0)
			{
				points[label - 1].push_back(Point(col, row));
			}
		}
	}

	// 3 检测结果
	vector<Result> t_results;
	for (int i = 0; i < connectNumber - 1; i++)
	{
		Result result;
		int label = i + 1;
		result.label = label;

		// 连通域面积 面积筛选
		int area = points[i].size();
		if (m_params.filterByArea)
		{
			if (area < m_params.minArea || area > m_params.maxArea)
			{
				continue;
			}
		}
		result.area = area;

		// 最小外接矩形 轴比筛选 矩形度筛选
		cv::RotatedRect minRect = cv::minAreaRect(points[i]);
		result.minRect = minRect;
		float axleratio = (float)min(minRect.size.height, minRect.size.width) / max(minRect.size.height, minRect.size.width);
		if (m_params.filterByInertia)
		{
			if (axleratio < m_params.minInertiaRatio || axleratio > m_params.maxInertiaRatio)
			{
				continue;
			}
		}
		result.axleratio = axleratio;
		

		// 外轮廓  计算连通域的周长 周长筛选 圆形度筛选
		cv::Mat componentMask = (m_connect == label);
		std::vector<std::vector<cv::Point>> contours;
		std::vector<cv::Point> externalContour;
		cv::findContours(componentMask, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

		for (const auto& contour : contours)
		{
			if (cv::contourArea(contour) > 0 && cv::contourArea(contour) < componentMask.total())
			{
				externalContour = contour;
			}
		}
		float perimeter = cv::arcLength(externalContour, true);
		if (m_params.filterByPerimeter)
		{
			if (perimeter < m_params.minPerimeter || perimeter > m_params.maxPerimeter)
			{
				continue;
			}
		}
		result.perimeter = perimeter;
		result.externalContour = externalContour;
		result.Contour = contours;

		Moments moms = moments(externalContour);
		float area2 = moms.m00;
		float circularity = 4 * CV_PI * area2 / (perimeter * perimeter);
		if (m_params.filterByCircularity)
		{
			if (circularity < m_params.minCircularity || circularity > m_params.maxCircularity)
				continue;
		}
		result.circularity = circularity;

		float rectangularity = (float)area2 / (minRect.size.height * minRect.size.width);
		if (m_params.filterByRectangularity)
		{
			if (rectangularity < m_params.minRectangularity || rectangularity > m_params.maxRectangularity)
			{
				continue;
			}
		}
		result.rectangularity = rectangularity;


		// 质心 质心偏移筛选
		float centroid_x = 0;
		float centroid_y = 0;
		for (auto& point : points[i])
		{
			centroid_x += point.x;
			centroid_y += point.y;
		}
		centroid_x /= points[i].size();
		centroid_y /= points[i].size();
		float center_x = minRect.center.x;
		float center_y = minRect.center.y;
		float centroidOffset = sqrt((centroid_x - center_x) * (centroid_x - center_x) + (centroid_y - center_y) * (centroid_y - center_y));
		if (m_params.filterByCentroidOffset)
		{
			if (centroidOffset < m_params.minCentroidOffset || circularity > m_params.maxCentroidOffset)
				continue;
		}
		result.centroid = Point2f(centroid_x, centroid_y);
		result.centroidOffset = centroidOffset;

		t_results.push_back(result);
	}

	// 4 排序
	if (m_params.sortMethod == 1)
	{
		switch (m_params.sortFeature)
		{
		case 0:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.area > b.area; });
			break;
		case 1:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.perimeter > b.perimeter; });
			break;
		case 2:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.circularity > b.circularity; });
			break;
		case 3:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.rectangularity > b.rectangularity; });
			break;
		case 4:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.axleratio > b.axleratio; });
			break;
		case 5:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.minRect.center.x > b.minRect.center.x; });
			break;
		case 6:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.minRect.center.y > b.minRect.center.y; });
			break;
		default:
			break;
		}
	}
	else if (m_params.sortMethod == 2)
	{
		switch (m_params.sortFeature)
		{
		case 0:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.area < b.area; });
			break;
		case 1:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.perimeter < b.perimeter; });
			break;
		case 2:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.circularity < b.circularity; });
			break;
		case 3:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.rectangularity < b.rectangularity; });
			break;
		case 4:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.axleratio < b.axleratio; });
			break;
		case 5:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.minRect.center.x < b.minRect.center.x; });
			break;
		case 6:
			std::sort(t_results.begin(), t_results.end(), [](Result a, Result b) {return a.minRect.center.y < b.minRect.center.y; });
			break;
		default:
			break;
		}
	}

	if (m_params.maxTargets)
	{
		if (t_results.size() > m_params.maxTargets)
		{
			t_results.resize(m_params.maxTargets);
		}
	}
	m_results = t_results;

	cv::Mat componentMask(m_connect.size(), CV_8UC1, Scalar(0));
	for (auto& result : m_results)
	{
		cv::Mat t_mask = (m_connect == result.label);
		cv::bitwise_or(componentMask, t_mask, componentMask);
	}
	m_resultBinary = componentMask.clone();

	// 结果显示
	if (m_show)
	{
		Mat show = Mat::zeros(m_binary.size(), CV_8UC3); //定义标记结果图像

		// 二值图像
		imshow("m_binary", m_binary);
		imshow("result binary", m_resultBinary);
		// 连通域图像
		{
			RNG rng(50); //RNG随机数生成器
			vector<Vec3b>colors; //定义颜色数组
			for (int i = 0; i < connectNumber; ++i) {
				//使用均匀分布的随机确定颜色
				Vec3b vec3 = Vec3b(rng.uniform(0, 256), rng.uniform(0, 256), rng.uniform(0, 256));
				colors.push_back(vec3);
			}

			//以不同颜色标记出不同的连通域
			int w = show.cols;
			int h = show.rows;
			for (int row = 0; row < h; ++row) {
				for (int col = 0; col < w; ++col) {
					int label = m_connect.at<uint32_t>(row, col);
					if (label == 0) {//背景的黑色不改变
						continue;
					}
					show.at<Vec3b>(row, col) = colors[label];
				}
			}
			//cout << "connectNumber: " << connectNumber << endl;
			Mat show0 = show.clone();
			cv::putText(show0, "connect number: " + to_string(connectNumber), Point(20, 20), 1, 1.5, Scalar(0, 0, 255), 2, LINE_8);
			cv::imshow("result", show0);
		}

		// 最小外接矩形图像
		{
			Mat show1 = show.clone();
			for (int j = 0; j < m_results.size(); j++)
			{
				cv::Point2f vertices[4];
				m_results[j].minRect.points(vertices);
				for (int i = 0; i < 4; i++)
				{
					cv::line(show1, vertices[i], vertices[(i + 1) % 4], cv::Scalar(0, 0, 255), 2);
				}
			}
			cv::putText(show1, "results number: " + to_string(m_results.size()), Point(20, 20), 1, 1.5, Scalar(0, 0, 255), 2, LINE_8);
			cv::imshow("minRect", show1);
		}


		// 所有轮廓图像
		{
			Mat show4 = show.clone();
			for (auto result : m_results)
			{
				circle(show4, result.minRect.center, 1, { 0, 255, 255 }, 2);
				circle(show4, result.centroid, 1, { 0, 0, 255 }, 2);
				for (auto Contour : result.Contour)
				{
					for (auto point : Contour)
					{
						circle(show4, point, 1, { 0, 0, 255 }, 1);
					}
				}
			}
			cv::putText(show4, "results number: " + to_string(m_results.size()), Point(20, 20), 1, 1.5, Scalar(0, 0, 255), 2, LINE_8);
			cv::imshow("contours", show4);
		}

		// 所有轮廓图像
		{
			Mat show3 = show.clone();
			for (auto result : m_results)
			{
				circle(show3, result.minRect.center, 1, { 0, 255, 255 }, 2);
				circle(show3, result.centroid, 1, { 0, 0, 255 }, 2);
				for (auto point : result.externalContour)
				{
					circle(show3, point, 1, { 0, 0, 255 }, 1);
				}
			}
			cv::putText(show3, "results number: " + to_string(m_results.size()), Point(20, 20), 1, 1.5, Scalar(0, 0, 255), 2, LINE_8);
			cv::imshow("externalContour", show3);
		}
	}
}


cv::Mat BlobDetector::getResults()
{
	if (m_results.empty())
	{
		return Mat(1, 1, CV_8UC1, Scalar(0));
	}
	Mat result = Mat(m_results.size(), 8, CV_32FC1, Scalar(0));

	for (int i = 0; i < m_results.size(); i++)
	{
		float* ptr = result.ptr<float>(i);
		ptr[0] = m_results[i].centroid.x;
		ptr[1] = m_results[i].centroid.y;
		ptr[2] = m_results[i].area;
		ptr[3] = m_results[i].perimeter;
		ptr[4] = m_results[i].circularity;
		ptr[5] = m_results[i].rectangularity;
		ptr[6] = m_results[i].axleratio;
		ptr[7] = m_results[i].centroidOffset;
	}

	return result;
}

cv::Mat BlobDetector::getMinRects()
{
	if (m_results.empty())
	{
		return Mat(1, 1, CV_8UC1, Scalar(0));
	}

	Mat result = Mat(m_results.size(), 13, CV_32FC1, Scalar(0));
	for (int i = 0; i < m_results.size(); i++)
	{
		float* ptr = result.ptr<float>(i);
		ptr[0] = m_results[i].minRect.center.x;
		ptr[1] = m_results[i].minRect.center.y;
		ptr[2] = m_results[i].minRect.angle;
		ptr[3] = max(m_results[i].minRect.size.height, m_results[i].minRect.size.width);
		ptr[4] = min(m_results[i].minRect.size.height, m_results[i].minRect.size.width);
		cv::Point2f vertices[4];
		m_results[i].minRect.points(vertices);
		ptr[5] = vertices[0].x;
		ptr[6] = vertices[0].y;
		ptr[7] = vertices[1].x;
		ptr[8] = vertices[1].y;
		ptr[9] = vertices[2].x;
		ptr[10] = vertices[2].y;
		ptr[11] = vertices[3].x;
		ptr[12] = vertices[3].y;
	}

	return result;
}

cv::Mat BlobDetector::getContours()
{
	if (m_results.empty())
	{
		return Mat(1, 1, CV_8UC1, Scalar(0));
	}
	int maxnum = 0;
	std::vector<int> pnum(m_results.size());
	for (int i = 0; i < m_results.size(); i++)
	{
		int num = 0;
		for (int j = 0; j < m_results[i].Contour.size(); j++)
		{
			num += m_results[i].Contour[j].size();
		}
		pnum[i] = num;
		maxnum = max(maxnum, num);
	}

	Mat result = Mat(m_results.size(), maxnum * 2 + 1, CV_32FC1, Scalar(0));
	for (int i = 0; i < m_results.size(); i++)
	{
		float* ptr = result.ptr<float>(i);
		ptr[0] = pnum[i];
		std::vector<std::vector<cv::Point>> Contour = m_results[i].Contour;
		int k = 0;
		for (int m = 0; m < Contour.size(); m++)
		{
			std::vector<cv::Point> tContour = Contour[m];
			for (int n = 0; n < tContour.size(); n++, k++)
			{
				ptr[2 * k + 1] = tContour[n].x;
				ptr[2 * k + 2] = tContour[n].y;
			}
		}
	}
	return result;
}

cv::Mat BlobDetector::getExternalContours()
{
	if (m_results.empty())
	{
		return Mat(1, 1, CV_8UC1, Scalar(0));
	}
	int maxnum = 0;
	for (int i = 0; i < m_results.size(); i++)
	{
		maxnum = max(maxnum, (int)m_results[i].externalContour.size());
	}
	Mat result = Mat(m_results.size(), maxnum * 2 + 1, CV_32FC1, Scalar(0));
	for (int i = 0; i < m_results.size(); i++)
	{
		float* ptr = result.ptr<float>(i);
		std::vector<cv::Point> externalContour = m_results[i].externalContour;
		ptr[0] = externalContour.size();
		for (int j = 0; j < externalContour.size(); j++)
		{
			ptr[2 * j + 1] = externalContour[j].x;
			ptr[2 * j + 2] = externalContour[j].y;
		}
	}
	return result;
}



