#include"CircleDetect.h"
#include"EDPF.h"
#include<Eigen/Dense>

using namespace Eigen;
using namespace cv;
using namespace std;


void CircleDetect::setSplit(float lineLen, float sharpAngle, float SegmentCurvature)
{
	m_params.lineLen = lineLen;
	m_params.sharpAngle = sharpAngle;
	m_params.SegmentCurvature = SegmentCurvature;
}

void CircleDetect::setInlier(float inlier, float closedInlier)
{
	m_params.inlier = inlier;
	m_params.closedInlier = closedInlier;
}

void CircleDetect::setDistance(float centerDis, float radiusDis)
{
	m_params.centerDis = centerDis;
	m_params.radiusDis = radiusDis;
}

void CircleDetect::setFiltSize(int filtSize)
{
	m_params.filtSize = filtSize;
}

void CircleDetect::setMinLen(float minLen)
{
	m_params.minLen = minLen;
}


// 2 检测闭合边缘段
void CircleDetect::extractClosedEdges(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& out_sigments)
{
	out_sigments.clear();
	Point st, ed;
	double dist;
	for (int i = 0; i < (in_sigments).size(); i++)
	{
		st = in_sigments[i].front();
		ed = in_sigments[i].back();
		dist = sqrt((st.x - ed.x) * (st.x - ed.x) + (st.y - ed.y) * (st.y - ed.y));
		if (dist <= 5)
		{
			out_sigments.push_back(in_sigments[i]);
		}
	}
}

// 3 以直代曲
double PerpendicularDistance(const Point& pt, const Point& lineStart, const Point& lineEnd)
{
	double dx = (double)lineEnd.x - lineStart.x;
	double dy = (double)lineEnd.y - lineStart.y;
	double mag = sqrt(dx * dx + dy * dy);
	if (mag > 0.0)
	{
		dx /= mag; dy /= mag;
	}
	double pvx = (double)pt.x - lineStart.x;
	double pvy = (double)pt.y - lineStart.y;
	double pvdot = dx * pvx + dy * pvy;
	double ax = pvx - pvdot * dx;
	double ay = pvy - pvdot * dy;
	return sqrt(ax * ax + ay * ay);
}

void RamerDouglasPeucker(const vector<Point>& pointList, double epsilon, vector<Point>& out)
{
	double dmax = 0.0;
	int index = 0;
	int end = pointList.size() - 1;
	for (int i = 1; i < end; i++)
	{
		double d = PerpendicularDistance(pointList[i], pointList[0], pointList[end]);
		if (d > dmax)
		{
			index = i;
			dmax = d;
		}
	}

	if (dmax > epsilon)
	{
		vector<Point> recResults1;
		vector<Point> recResults2;
		vector<Point> firstLine(pointList.begin(), pointList.begin() + index + 1);
		vector<Point> lastLine(pointList.begin() + index, pointList.end());
		RamerDouglasPeucker(firstLine, epsilon, recResults1);
		RamerDouglasPeucker(lastLine, epsilon, recResults2);

		// Build the result list
		out.assign(recResults1.begin(), recResults1.end() - 1);
		out.insert(out.end(), recResults2.begin(), recResults2.end());
		if (out.size() < 2)
			throw runtime_error("Problem assembling output");
	}
	else
	{
		//Just return start and end points
		out.clear();
		out.push_back(pointList[0]);//pointList[0]
		out.push_back(pointList[end]);
	}
}

void CircleDetect::RDP(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& out_sigments)
{
	out_sigments.clear();
	for (auto& seg : in_sigments)
	{
		vector<Point> segTemp;

		double dmax = 0.0;
		int index = 0;
		int end = seg.size() - 1;

		for (int i = 1; i < end; i++)
		{
			double d = PerpendicularDistance(seg[i], seg[0], seg[end]);
			if (d > dmax)
			{
				index = i;
				dmax = d;
			}
		}

		// If max distance is greater than epsilon, recursively simplify
		float epsilon = 3;
		if (dmax > epsilon)
		{
			segTemp.clear();
			// Recursive call
			vector<Point> recResults1;
			vector<Point> recResults2;
			vector<Point> firstLine(seg.begin(), seg.begin() + index + 1);
			vector<Point> lastLine(seg.begin() + index, seg.end());
			RamerDouglasPeucker(firstLine, epsilon, recResults1);
			RamerDouglasPeucker(lastLine, epsilon, recResults2);
			
			// Build the result list
			segTemp.assign(recResults1.begin(), recResults1.end() - 1);
			segTemp.insert(segTemp.end(), recResults2.begin(), recResults2.end());
			if (segTemp.size() < 2)
				throw runtime_error("Problem assembling output");;
		}
		else
		{
			
			//Just return start and end points
			segTemp.clear();
			segTemp.push_back(seg[0]);//pointList[0]
			segTemp.push_back(seg[end]);
		}
		out_sigments.push_back(segTemp);

	}
}

// 4 链中角度较大的拐点处断开
void CircleDetect::rejectSharpTurn(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
	std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists)
{
	int no_seg_grps = in_sigments.size();
	vector<Point>break_points;
	double Threshold_theta = CV_PI / 2;
	for (int ii = 0; ii < no_seg_grps; ii++)//no_seg_grps
	{

		vector<Point> present_seg_grp = in_sigments[ii];
		int no_of_seg = present_seg_grp.size() - 1;//seg_point -1=no_of_seg;
		// use cos(angle) calculation to determine sharp turn angles
		Point present_vector = Point(present_seg_grp[1].x - present_seg_grp[0].x, present_seg_grp[1].y - present_seg_grp[0].y);

		for (int seg_no = 0; seg_no < no_of_seg - 1; seg_no++)
		{
			double length_present_vector = sqrt(present_vector.x * present_vector.x + present_vector.y * present_vector.y);//两点距离
			Point next_vector = Point(present_seg_grp[seg_no + 2].x - present_seg_grp[seg_no + 1].x, present_seg_grp[seg_no + 2].y - present_seg_grp[seg_no + 1].y);
			double length_next_vector = sqrt(next_vector.x * next_vector.x + next_vector.y * next_vector.y);
			//余弦定理
			double cos_pre_next = ((double)present_vector.x * next_vector.x + (double)present_vector.y * next_vector.y) / (length_present_vector * length_next_vector);

			if (cos_pre_next <= cos(m_params.sharpAngle / 180.0 * CV_PI))
			{
				break_points.push_back(Point(ii, seg_no + 1));//记录拐点(链编号，第几个点)
			}
			present_vector = next_vector;
		}//endfor
	}
	if (break_points.empty())//no break points
	{
		out_edgelists = in_edgelists;
		out_sigments = in_sigments;
		return;
	}
	int index = 0;
	int current_break = break_points[index].x;//链编号
	for (int ii = 0; ii < no_seg_grps; ii++)
	{
		vector<Point> current_seg = in_sigments[ii];
		vector<Point> current_edge = in_edgelists[ii];
		if (current_seg.size() > 2)
		{
			if (ii == current_break)
			{
				int count = 1;
				int first_edge_index, last_edge_index, first_seg_index, last_seg_index;
				while (ii == current_break)
				{
					if (count == 1)
					{
						first_seg_index = 0;
						first_edge_index = 0;
					}
					else
					{
						first_seg_index = last_seg_index;
						first_edge_index = last_edge_index;
					}
					last_seg_index = break_points[index].y;//链中拐点的编号

					for (int jj = first_edge_index; jj < current_edge.size(); jj++)
					{
						//在current_edge中找到拐点
						if ((current_seg[last_seg_index].x == current_edge[jj].x) && (current_seg[last_seg_index].y == current_edge[jj].y))
						{
							last_edge_index = jj;
							break;
						}//endif
					}//endfor

					//拐点前的链部分
					if ((last_seg_index - first_seg_index) >= 1)//连续两个拐点的编号距离
					{
						vector<Point> block_seg, block_edge;
						for (int kk = first_seg_index; kk <= last_seg_index; kk++)
						{
							block_seg.push_back(current_seg[kk]);
						}
						for (int kk = first_edge_index; kk <= last_edge_index; kk++)
						{
							block_edge.push_back(current_edge[kk]);
						}
						out_edgelists.push_back(block_edge);
						out_sigments.push_back(block_seg);
					}

					//check for break
					index += 1;
					if (index > (break_points.size() - 1))
					{
						break;
					}
					current_break = break_points[index].x;
					count = count + 1;

				}//endwhile

				//拐点后的链部分
				if ((current_seg.size() - last_seg_index) >= 2)
				{
					vector<Point> block1_seg, block1_edge;
					for (int ll = last_seg_index; ll < current_seg.size(); ll++)
					{

						block1_seg.push_back(current_seg[ll]);

					}
					for (int ll = last_edge_index; ll < current_edge.size(); ll++)
					{

						block1_edge.push_back(current_edge[ll]);
					}
					out_edgelists.push_back(block1_edge);
					out_sigments.push_back(block1_seg);

				}
			}//endif
			else
			{
				out_edgelists.push_back(in_edgelists[ii]);
				out_sigments.push_back(in_sigments[ii]);

			}//endelse
		}//endif
	}//endfor
}

// 5 直线弧线交界处断开
void CircleDetect::detectline(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
							  std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists)
{
	int no_seg_grps = in_sigments.size();
	vector<Point> break_points;

	for (int ii = 0; ii < no_seg_grps; ii++)//no_seg_grps
	{
		vector<Point> present_seg_grp = in_sigments[ii];
		int no_of_seg = present_seg_grp.size();//seg_point -1=no_of_seg;

		for (int seg_no = 0; seg_no < no_of_seg - 1; seg_no++)
		{
			double length = sqrt((present_seg_grp[seg_no + 1].x - present_seg_grp[seg_no].x) * (present_seg_grp[seg_no + 1].x - present_seg_grp[seg_no].x)
				+ (present_seg_grp[seg_no + 1].y - present_seg_grp[seg_no].y) * (present_seg_grp[seg_no + 1].y - present_seg_grp[seg_no].y));//两点距离
			if ((int)length > m_params.lineLen)
			{
				if (seg_no != 0)
				{
					Point temp = Point(ii, seg_no);
					break_points.push_back(temp);
				}
				if (seg_no + 1 != no_of_seg - 1)
				{
					Point temp = Point(ii, seg_no + 1);
					break_points.push_back(temp);
				}
			}

		}
	}

	if (break_points.empty())//no break points
	{
		out_edgelists = in_edgelists;
		out_sigments = in_sigments;
	}
	else
	{
		int index = 0;
		int current_break = break_points[index].x;//链编号
		for (int ii = 0; ii < no_seg_grps; ii++)
		{
			vector<Point> current_seg = in_sigments[ii];
			vector<Point> current_edge = in_edgelists[ii];
			if (ii == current_break)
			{
				int count = 1;
				int first_edge_index, last_edge_index, first_seg_index, last_seg_index;
				while (ii == current_break)
				{
					if (count == 1)
					{
						first_seg_index = 0;
						first_edge_index = 0;
					}
					else
					{
						first_seg_index = last_seg_index;
						first_edge_index = last_edge_index;
					}
					last_seg_index = break_points[index].y;//链中拐点的编号

					for (int jj = first_edge_index; jj < current_edge.size(); jj++)
					{
						//在current_edge中找到拐点
						if ((current_seg[last_seg_index].x == current_edge[jj].x) && (current_seg[last_seg_index].y == current_edge[jj].y))
						{
							last_edge_index = jj;
							break;
						}//endif
					}//endfor
					//block before break
					//cut block



					//拐点前的链部分
					if ((last_seg_index - first_seg_index) >= 1)//连续两个拐点的编号距离
					{
						vector<Point> block_seg, block_edge;
						for (int kk = first_seg_index; kk <= last_seg_index; kk++)
						{
							block_seg.push_back(current_seg[kk]);
						}
						for (int kk = first_edge_index; kk <= last_edge_index; kk++)
						{
							block_edge.push_back(current_edge[kk]);
						}
						out_edgelists.push_back(block_edge);
						out_sigments.push_back(block_seg);
					}

					//check for break
					index += 1;
					if (index > (break_points.size() - 1))
					{
						break;
					}
					current_break = break_points[index].x;
					count = count + 1;

				}//endwhile
				//block after break
				//拐点后的链部分
				if ((current_seg.size() - last_seg_index) >= 2)
				{
					vector<Point> block1_seg, block1_edge;
					for (int ll = last_seg_index; ll < current_seg.size(); ll++)
					{

						block1_seg.push_back(current_seg[ll]);

					}
					for (int ll = last_edge_index; ll < current_edge.size(); ll++)
					{

						block1_edge.push_back(current_edge[ll]);
					}
					out_edgelists.push_back(block1_edge);
					out_sigments.push_back(block1_seg);

				}
			}//endif
			else
			{
				out_edgelists.push_back(in_edgelists[ii]);
				out_sigments.push_back(in_sigments[ii]);

			}//endelse
		}//endfor
	}
}

// 6 根据链的走向、弯曲方向变化断开链
void CircleDetect::detectInflexPt(std::vector<std::vector<cv::Point>>& in_sigments, std::vector<std::vector<cv::Point>>& in_edgelists,
								  std::vector<std::vector<cv::Point>>& out_sigments, std::vector<std::vector<cv::Point>>& out_edgelists)
{
	int no_seg_grps = in_sigments.size();
	for (int ii = 0; ii < no_seg_grps; ii++)
	{
		vector<Point> present_seg_grp;
		vector<int> break_at_seg, break_points_y_x_list;
		present_seg_grp = in_sigments[ii];
		if (present_seg_grp.size() <= 4)
		{

			out_sigments.push_back(in_sigments[ii]);
			out_edgelists.push_back(in_edgelists[ii]);
			continue;
		}//endif

		// slope angle	
		vector<float> theta, theta_2pi;
		for (int jj = 0; jj < present_seg_grp.size() - 1; jj++)// adjoint angles
		{
			//角度
			float tempAngle = atan2(present_seg_grp[jj + 1].x - present_seg_grp[jj].x, present_seg_grp[jj + 1].y - present_seg_grp[jj].y);
			theta.push_back(tempAngle);
			if (tempAngle < 0)//顺时针
			{
				theta_2pi.push_back(tempAngle + 2 * CV_PI);//转换为正数
			}//endif
			else//逆时针
			{
				theta_2pi.push_back(tempAngle);
			}
		}//endfor

		// angle differences
		vector<float> theta_diff, theta_diff_2pi;
		for (int kk = 0; kk < theta.size() - 1; kk++)//两个角度之差
		{
			theta_diff.push_back(theta[kk] - theta[kk + 1]);
			theta_diff_2pi.push_back(theta_2pi[kk] - theta_2pi[kk + 1]);
		}//endfor
		// correct angles
		for (int i = 0; i < theta_diff.size(); i++)//两个角度之差
		{
			if (theta_diff[i]<-1 * CV_PI || theta_diff[i]> CV_PI)
			{
				theta_diff[i] = theta_diff_2pi[i];
			}//endif
		}//endfor
		//polarity
		vector<int> polarity;
		for (int i = 0; i < theta_diff.size(); i++)
		{
			if (theta_diff[i] > 0)
			{
				polarity.push_back(1);
			}//endif
			else { polarity.push_back(0); }//endelse
		}//endfor


		int count = 0;// non zero polarity
		int polaritySize = polarity.size();
		for (int i = 0; i < polaritySize; i++)
		{
			if (polarity[i] == 1) count++;
		}//endfor

		if (count > (double)(polaritySize / 2))
		{

			for (int j = 0; j < polaritySize; j++)
			{
				if (polarity[j] == 0)
				{

					polarity[j] = 1;
				}
				else
				{
					polarity[j] = 0;
				}
			}//endfor
		}//endif


		//checking ...0 1 0...type inflexion
		vector<int> mask1;
		mask1.push_back(0);
		mask1.push_back(1);
		mask1.push_back(0);
		int location = 0;
		while (location < (polaritySize - 2))
		{

			vector<int> window;
			window.push_back(polarity[location]);
			window.push_back(polarity[location + 1]);
			window.push_back(polarity[location + 2]);
			vector<int> is_inflexion;
			// xor, result is stored in is_inflexion
			for (int i = 0; i < 3; i++)
			{
				if (window[i] - mask1[i] == 0) { is_inflexion.push_back(0); }
				else
				{
					is_inflexion.push_back(1);
				}
			}//endfor
			vector<int>::iterator pos;
			pos = find(is_inflexion.begin(), is_inflexion.end(), 1);
			if (pos != is_inflexion.end())
			{
				location += 1;
			}
			else
			{
				polarity[location] = 0;
				polarity[location + 1] = 0;
				polarity[location + 2] = 0;
				break_at_seg.push_back(location + 2);
				location += 2;
			}//endifelse
		}//endwhile


		//checking ... 0 0 1 1 ...type inflexion
		vector<int> mask2;
		mask2.push_back(0);
		mask2.push_back(1);
		mask2.push_back(1);
		mask2.push_back(0);
		int location2 = 0;
		while (location2 < polaritySize - 3)
		{
			vector<int> window;
			window.push_back(polarity[location2]);
			window.push_back(polarity[location2 + 1]);
			window.push_back(polarity[location2 + 2]);
			window.push_back(polarity[location2 + 3]);
			vector<int> is_inflexion;
			// xor, result is stored in is_inflexion
			for (int i = 0; i < 4; i++)
			{
				if (window[i] - mask2[i] == 0) { is_inflexion.push_back(0); }
				else {
					is_inflexion.push_back(1);
				}
			}//endfor
			vector<int>::iterator pos;
			pos = find(is_inflexion.begin(), is_inflexion.end(), 1);
			if (pos != is_inflexion.end())
			{
				location2 += 1;
			}
			else
			{

				polarity[location2] = 0;
				polarity[location2 + 1] = 0;
				polarity[location2 + 2] = 0;
				polarity[location2 + 3] = 0;
				break_at_seg.push_back(location2 + 2);
				break_at_seg.push_back(location2 + 3);
				location2 += 3;
			}//endifelse
		}//endwhile

		vector<int> transitions;
		vector<int> locations;
		// xor, result is stored in is_inflexion
		for (int i = 0; i < polaritySize - 1; i++)
		{
			if ((polarity[i] - polarity[i + 1]) == 0)
			{
				transitions.push_back(0);
			}
			else {
				transitions.push_back(1);
				locations.push_back(i); // record the index of nonzero elements
			}
		}//endfor

		// checking 1 0... type inflexion (the start segment is off)
		if ((!locations.empty()) && (locations.front() == 0))
		{

			break_at_seg.push_back(1);
			locations.erase(locations.begin());

		}//endif

		//checking ... 0 1 type inflexion(the end segment is off)
		if ((!locations.empty()) && (locations.back() == (polaritySize - 2)))//*(locations.end()-1
		{

			break_at_seg.push_back(present_seg_grp.size() - 2);
			locations.pop_back();
		}

		//checking ...0 1 1 1...type inflexion
		if (!locations.empty())
		{
			for (auto it = locations.begin(); it != locations.end(); it++)
				break_at_seg.push_back((*it) + 1);
		}//endif

		//now_breaking the edges and seglists
		if (break_at_seg.empty())
		{
			out_sigments.push_back(in_sigments[ii]);
			out_edgelists.push_back(in_edgelists[ii]);
		}//endif
		else
		{
			sort(break_at_seg.begin(), break_at_seg.end());// sort in the ascending order
			int breakAtSegSize = break_at_seg.size();
			for (int jj = 0; jj < breakAtSegSize; jj++)
			{

				int x1 = present_seg_grp[break_at_seg[jj]].x;
				int y1 = present_seg_grp[break_at_seg[jj]].y;

				vector<int> a;

				for (int kk = 0; kk < in_edgelists[ii].size(); kk++)
				{
					if ((in_edgelists[ii][kk].x == x1) && (in_edgelists[ii][kk].y == y1))
					{

						a.push_back(kk);
					}//endif
				}//endfor
				if (!a.empty())
				{

					break_points_y_x_list.push_back(a[0]);
				}//endif
			}//endfor



			// for segList
			vector<vector<Point> >seglist_temp;
			int k = 0;
			for (int jj = 0; jj < breakAtSegSize; jj++)
			{
				vector<Point> seg_temp;
				for (int p = k; p <= break_at_seg[jj]; p++)
				{
					seg_temp.push_back(in_sigments[ii][p]);
				}
				seglist_temp.push_back(seg_temp);
				k = break_at_seg[jj];
			}
			vector<Point> seg_temp;
			for (int pp = k; pp < in_sigments[ii].size(); pp++)// the final segment
			{
				seg_temp.push_back(in_sigments[ii][pp]);
			}//endfor
			seglist_temp.push_back(seg_temp);
			for (int ppp = 0; ppp < seglist_temp.size(); ppp++)
			{
				out_sigments.push_back(seglist_temp[ppp]);
			}//endfor

			// for edgeList
			vector<vector<Point> >edgelist_temp;
			int k2 = 0;
			if (!break_points_y_x_list.empty())
			{
				for (int jj = 0; jj < break_points_y_x_list.size(); jj++)
				{
					vector<Point> edge_temp;
					for (int p = k2; p <= break_points_y_x_list[jj]; p++)
					{
						edge_temp.push_back(in_edgelists[ii][p]);
					}
					edgelist_temp.push_back(edge_temp);
					k2 = break_points_y_x_list[jj];
				}
				vector<Point> edge_temp2;
				for (int pp = k2; pp < in_edgelists[ii].size(); pp++)// the final segment
				{
					edge_temp2.push_back(in_edgelists[ii][pp]);
				}//endfor
				edgelist_temp.push_back(edge_temp2);
				for (int ppp = 0; ppp < edgelist_temp.size(); ppp++)
				{
					out_edgelists.push_back(edgelist_temp[ppp]);
				}//endfor
			}//endif
		}//endelse
	}
}


// 7 删除段的链和直线
void CircleDetect::detectshort(std::vector<std::vector<cv::Point>>& in_out_edgelists)
{
	vector<vector<Point>>::iterator it = in_out_edgelists.begin();
	while (it != in_out_edgelists.end())
	{
		/*compute the line segment generated by the two endpoints of the arc,
		and then judge the midpoint of the arc if lying on or near the line
		*/
		Point edgeSt = Point((*it).front().x, (*it).front().y);
		Point edgeEd = Point((*it).back().x, (*it).back().y);
		int midIndex = (*it).size() / 2;

		Point edgeMid = Point((*it)[midIndex].x, (*it)[midIndex].y);

		double distStEd = sqrt((edgeSt.x - edgeEd.x) * (edgeSt.x - edgeEd.x) + (edgeSt.y - edgeEd.y) * (edgeSt.y - edgeEd.y));
		double distStMid = sqrt((edgeSt.x - edgeMid.x) * (edgeSt.x - edgeMid.x) + (edgeSt.y - edgeMid.y) * (edgeSt.y - edgeMid.y));
		double distMidEd = sqrt((edgeEd.x - edgeMid.x) * (edgeEd.x - edgeMid.x) + (edgeEd.y - edgeMid.y) * (edgeEd.y - edgeMid.y));
		double distDifference = abs((distStMid + distMidEd) - distStEd);


		if ((*it).size() <= m_params.minLen || distDifference <= m_params.SegmentCurvature * (distStMid + distMidEd))// 2 3 fixed number; (*it).size() <=20
		{
			it = in_out_edgelists.erase(it);
		}
		else { it++; }
	}//endwhile
}

// 8 提取闭合边和非闭合边
void CircleDetect::extractClosedEdges2(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& out_edgelists
						, std::vector<std::vector<cv::Point>>& out_closededgelists)
{
	Point st, ed;
	double dist;
	for (int i = 0; i < (in_edgelists).size(); i++)
	{
		st = in_edgelists[i].front();
		ed = in_edgelists[i].back();
		dist = sqrt((st.x - ed.x) * (st.x - ed.x) + (st.y - ed.y) * (st.y - ed.y));
		if (dist <= 5)
		{
			out_closededgelists.push_back(in_edgelists[i]);
		}
		else
		{
			out_edgelists.push_back(in_edgelists[i]);
		}
	}
}

//8 非闭合曲线分组
inline int Sign(float x)
{
	int sign;
	if (x > 0)
		return 1;
	else
		return -1;
}
bool twoArcsCenterRadius(vector<Point>& A1B1C1, vector<Point>& A2B2C2, bool* flag, Vec3f* temp1_center_radius, Vec3f* temp2_center_radius, float T_o, float T_r)
{
	*flag = false;
	// constraint 1:
	//S2 E2 and M1 lie in different sizes of Line S1E1 && S1 E1 and M2 lie in different sizes of Line S2E2
	Point S1 = A1B1C1[2];//A1B1C1.front()
	Point E1 = *(A1B1C1.end() - 2);//A1B1C1.back()
	Point S2 = A2B2C2[2];//A2B2C2.front()
	Point E2 = *(A2B2C2.end() - 2);//A2B2C2.back()
	Point M1 = A1B1C1[A1B1C1.size() / 2];
	Point M2 = A2B2C2[A2B2C2.size() / 2];
	// Line S1E1, S2E2
	//两段弧边界点的直线斜率
	float K_S1E1 = ((double)S1.x - E1.x) / ((double)S1.y - E1.y + 1e-6);
	float K_S2E2 = ((double)S2.x - E2.x) / ((double)S2.y - E2.y + 1e-6);
	// the location of midpoint M1 regarding line S1E1, M2 regarding line S2E2
	//弧中点与直线位置，Sign：取数字n的符号,大于0返回1,小于0返回-1,等于0返回0，在直线的上方还是下方
	int SignM1 = Sign(M1.x - S1.x - K_S1E1 * (M1.y - S1.y));
	int SignM2 = Sign(M2.x - S2.x - K_S2E2 * (M2.y - S2.y));
	// S2, E2 of A2B2C2 should lie in different sides of M1 && S1, E1 of A1B1C1 should lie in different sides of M2
	//直线应夹在中点和另一段弧之间
	int SignS2 = Sign(S2.x - S1.x - K_S1E1 * (S2.y - S1.y));//s2相对直线1位置
	int SignE2 = Sign(E2.x - S1.x - K_S1E1 * (E2.y - S1.y));//e2相对直线1位置
	int SignS1 = Sign(S1.x - S2.x - K_S2E2 * (S1.y - S2.y));//s1相对直线2位置
	int SignE1 = Sign(E1.x - S2.x - K_S2E2 * (E1.y - S2.y));//e1相对直线2位置

	//如果起点终点在同一侧
	if (SignS1 * SignE1 >= 0 || SignS2 * SignE2 >= 0)// different sides then check constraint 2
	{
		//if constratint 1 holds, then check constratint 2
		// 
		//去除两段弧前后5个点
		Point A1 = A1B1C1[3];//A1B1C1.front()  A1B1C1[5]
		Point C1 = *(A1B1C1.end() - 4);//A1B1C1.back() *(A1B1C1.end()-6)
		int firstMidIndex = A1B1C1.size() / 2;
		Point B1 = A1B1C1[firstMidIndex];

		Point A2 = *(A2B2C2.begin() + 3);//A2B2C2.front()
		Point C2 = *(A2B2C2.end() - 4);//A2B2C2.back()
		int secondMidIndex = A2B2C2.size() / 2;
		Point B2 = A2B2C2[secondMidIndex];


		// 裁剪后的弧的弦长度
		double A1C1 = sqrt((A1.x - C1.x) * (A1.x - C1.x) + (A1.y - C1.y) * (A1.y - C1.y));
		double A2C2 = sqrt((A2.x - C2.x) * (A2.x - C2.x) + (A2.y - C2.y) * (A2.y - C2.y));
		// use the co-circle theorem to verify if the first arc and the second arc are 
		// from the same circle
		//the length for triangle edges
		//两个三角形的边长度
		double A1B2 = sqrt((A1.x - B2.x) * (A1.x - B2.x) + (A1.y - B2.y) * (A1.y - B2.y));
		double C1B2 = sqrt((C1.x - B2.x) * (C1.x - B2.x) + (C1.y - B2.y) * (C1.y - B2.y));
		double A2B1 = sqrt((A2.x - B1.x) * (A2.x - B1.x) + (A2.y - B1.y) * (A2.y - B1.y));
		double C2B1 = sqrt((C2.x - B1.x) * (C2.x - B1.x) + (C2.y - B1.y) * (C2.y - B1.y));

		// the areas of the two triangles两个三角形的面积  海伦公式sqrt(p*(p-a)*(p-b)*(p-c)),p为周长的一半
		double p1 = (A1B2 + C1B2 + A1C1) / 2;//周长的一半
		double S_A1B2C1 = sqrt(p1 * (p1 - A1B2) * (p1 - C1B2) * (p1 - A1C1));// Hallen formulation 海伦公式求三角形面积
		double p2 = (A2B1 + C2B1 + A2C2) / 2;
		double S_A2B1C2 = sqrt(p2 * (p2 - A2B1) * (p2 - C2B1) * (p2 - A2C2));// Hallen formulation 

		// use the ratio of the multiplication of edge length to area
		//三角形外界圆的半径：r=abc/4s，公式中a，b，c分别为三角形的三边，s为面积。
		double R1 = (A1B2 * C1B2 * A1C1) / (4 * S_A1B2C1);
		double R2 = (A2B1 * C2B1 * A2C2) / (4 * S_A2B1C2);

		if (abs(R1 - R2) <= T_r)//less than some certain threshold->group
		{
			//calculate the circumcircle center from the triangles' vertices 
			//根据三角形求圆心，圆心与三个点的距离等于半径
			// for A1B2C1
			double a11 = 2 * ((double)B2.y - A1.y);
			double b11 = 2 * ((double)B2.x - A1.x);
			double c11 = (double)B2.y * B2.y + B2.x * B2.x - A1.y * A1.y - A1.x * A1.x;
			double a12 = 2 * ((double)C1.y - B2.y);
			double b12 = 2 * ((double)C1.x - B2.x);
			double c12 = (double)C1.y * C1.y + C1.x * C1.x - B2.y * B2.y - B2.x * B2.x;

			//for A2B1C2
			double a21 = (double)2 * (B1.y - A2.y);
			double b21 = (double)2 * (B1.x - A2.x);
			double c21 = B1.y * B1.y + B1.x * B1.x - A2.y * A2.y - A2.x * A2.x;
			double a22 = (double)2 * (C2.y - B1.y);
			double b22 = 2 * (C2.x - B1.x);
			double c22 = C2.y * C2.y + C2.x * C2.x - B1.y * B1.y - B1.x * B1.x;

			//centers求圆心
			Point2f O1, O2;
			O1.x = (a11 * c12 - a12 * c11) / (a11 * b12 - a12 * b11);
			O1.y = (c11 * b12 - c12 * b11) / (a11 * b12 - a12 * b11);

			O2.x = (a21 * c22 - a22 * c21) / (a21 * b22 - a22 * b21);
			O2.y = (c21 * b22 - c22 * b21) / (a21 * b22 - a22 * b21);

			//两个圆心距离
			double distO1O2 = sqrt((O1.y - O2.y) * (O1.y - O2.y) + (O1.x - O2.x) * (O1.x - O2.x));

			if (distO1O2 <= T_o)
			{

				// compute inlier ratio
				double R = (R1 + R2) / 2;//圆半径为两个弧半径的均值
				Point2f O;
				//圆心：
				O.x = (O1.x + O2.x) / 2;
				O.y = (O1.y + O2.y) / 2;

				int num1 = 0, num2 = 0;
				double dis1, dis2;
				for (auto it1 = A1B1C1.begin(); it1 != A1B1C1.end(); it1++)
				{
					dis1 = sqrt(((*it1).x - O.x) * ((*it1).x - O.x) + ((*it1).y - O.y) * ((*it1).y - O.y));//弧上的点与圆心的距离
					if (abs(dis1 - R) <= 2)
						num1++;//弧1在圆上点的数量
				}//endfor

				for (auto it2 = A2B2C2.begin(); it2 != A2B2C2.end(); it2++)
				{
					dis2 = sqrt(((*it2).x - O.x) * ((*it2).x - O.x) + ((*it2).y - O.y) * ((*it2).y - O.y));
					if (abs(dis2 - R) <= 2)
						num2++;//弧2在圆上点的数量
				}//endfor
				int size1 = A1B1C1.size();
				int size2 = A2B2C2.size();
				double edgeInlier = (double)(num1 + num2) / (size1 + size2);//弧在圆上点的数量与弧总点数的比率

				if (edgeInlier >= 0.2)//0.4
				{
					*flag = true;
					(*temp1_center_radius)[0] = O1.x;
					(*temp1_center_radius)[1] = O1.y;
					(*temp1_center_radius)[2] = R1;
					(*temp2_center_radius)[0] = O2.x;
					(*temp2_center_radius)[1] = O2.y;
					(*temp2_center_radius)[2] = R2;
				}
			}//endif
		}
	}
	return true;
}

bool comCirCenterRadius(Point& A, Point& B, Point& C, double* R, Point2f* O)
{

	//the length for triangle edges
	double AB = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
	double CB = sqrt((C.x - B.x) * (C.x - B.x) + (C.y - B.y) * (C.y - B.y));
	double AC = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y));
	// the areas of the two triangles
	double p = (AB + CB + AC) / 2;
	double S_ABC = sqrt(p * (p - AB) * (p - CB) * (p - AC));// Hallen formulation 
	// radius from areas
	*R = (AB * CB * AC) / (4 * S_ABC);

	//calculate the circular center based on the triangles' vertices
	// for ABC
	double a11 = 2 * (B.y - A.y);
	double b11 = 2 * (B.x - A.x);
	double c11 = B.y * B.y + B.x * B.x - A.y * A.y - A.x * A.x;
	double a12 = 2 * (C.y - B.y);
	double b12 = 2 * (C.x - B.x);
	double c12 = C.y * C.y + C.x * C.x - B.y * B.y - B.x * B.x;
	(*O).x = (a11 * c12 - a12 * c11) / (a11 * b12 - a12 * b11);
	(*O).y = (c11 * b12 - c12 * b11) / (a11 * b12 - a12 * b11);
	return true;

}
bool estimateSingleCenterRadius(vector<Point>& A1B1C1, double* estimateR, Point2f* estimateO)
{
	Point A1 = A1B1C1[5];//A1B1C1.front()
	Point C1 = *(A1B1C1.end() - 6);//A1B1C1.back()
	int firstMidIndex = A1B1C1.size() / 2;
	Point B1 = A1B1C1[firstMidIndex];
	Point D11 = A1B1C1[5 + (firstMidIndex - 5) / 2];
	Point D12 = A1B1C1[firstMidIndex + (A1B1C1.size() - firstMidIndex) / 2];


	// first estimate center and radius using arcs by themselves

	// for arc A1B1C1
	// S_A1B1C1
	double R_A1B1C1;
	Point2f O_A1B1C1;
	comCirCenterRadius(A1, B1, C1, &R_A1B1C1, &O_A1B1C1);

	// S_A1D11C1
	double R_A1D11C1;
	Point2f O_A1D11C1;
	comCirCenterRadius(A1, D11, C1, &R_A1D11C1, &O_A1D11C1);

	// S_A1D12C1
	double R_A1D12C1;
	Point2f O_A1D12C1;
	comCirCenterRadius(A1, D12, C1, &R_A1D12C1, &O_A1D12C1);


	// using median computes centers and radius
	vector<double> tempR, tempOX, tempOY;
	// R 
	tempR.push_back(R_A1B1C1);
	tempR.push_back(R_A1D11C1);
	tempR.push_back(R_A1D12C1);
	//OX
	tempOX.push_back(O_A1B1C1.x);
	tempOX.push_back(O_A1D11C1.x);
	tempOX.push_back(O_A1D12C1.x);
	//OY
	tempOY.push_back(O_A1B1C1.y);
	tempOY.push_back(O_A1D11C1.y);
	tempOY.push_back(O_A1D12C1.y);

	sort(tempR.begin(), tempR.end());
	*estimateR = tempR[1];

	sort(tempOX.begin(), tempOX.end());
	double estimateOX = tempOX[1];

	sort(tempOY.begin(), tempOY.end());
	double estimateOY = tempOY[1];
	*estimateO = Point2f(estimateOX, estimateOY);
	return true;

}

bool estimateCenterRadius(vector<Point>& A1B1C1, vector<Point>& A2B2C2, double* estimateR, Point2f* estimateO)
{
	Point A1 = A1B1C1[5];//A1B1C1.front()
	Point C1 = *(A1B1C1.end() - 6);//A1B1C1.back()
	int firstMidIndex = A1B1C1.size() / 2;
	Point B1 = A1B1C1[firstMidIndex];

	Point A2 = *(A2B2C2.begin() + 5);//A2B2C2.front()
	Point C2 = *(A2B2C2.end() - 6);//A2B2C2.back()
	int secondMidIndex = A2B2C2.size() / 2;
	Point B2 = A2B2C2[secondMidIndex];

	// the length of edge A1C1 
	double A1C1 = sqrt((A1.x - C1.x) * (A1.x - C1.x) + (A1.y - C1.y) * (A1.y - C1.y));
	double A2C2 = sqrt((A2.x - C2.x) * (A2.x - C2.x) + (A2.y - C2.y) * (A2.y - C2.y));

	//the length for triangle edges
	double A1B2 = sqrt((A1.x - B2.x) * (A1.x - B2.x) + (A1.y - B2.y) * (A1.y - B2.y));
	double C1B2 = sqrt((C1.x - B2.x) * (C1.x - B2.x) + (C1.y - B2.y) * (C1.y - B2.y));
	double A2B1 = sqrt((A2.x - B1.x) * (A2.x - B1.x) + (A2.y - B1.y) * (A2.y - B1.y));
	double C2B1 = sqrt((C2.x - B1.x) * (C2.x - B1.x) + (C2.y - B1.y) * (C2.y - B1.y));

	// the areas of the two triangles
	double p1 = (A1B2 + C1B2 + A1C1) / 2;
	double S_A1B2C1 = sqrt(p1 * (p1 - A1B2) * (p1 - C1B2) * (p1 - A1C1));// Hallen formulation 
	double p2 = (A2B1 + C2B1 + A2C2) / 2;
	double S_A2B1C2 = sqrt(p2 * (p2 - A2B1) * (p2 - C2B1) * (p2 - A2C2));// Hallen formulation 

	// use the ratio of the multiplication of edge length to area
	double R1 = (A1B2 * C1B2 * A1C1) / (4 * S_A1B2C1);
	double R2 = (A2B1 * C2B1 * A2C2) / (4 * S_A2B1C2);

	//calculate the circumcircle center from the triangles' vertices
	// for A1B2C1
	double a11 = 2 * (B2.y - A1.y);
	double b11 = 2 * (B2.x - A1.x);
	double c11 = B2.y * B2.y + B2.x * B2.x - A1.y * A1.y - A1.x * A1.x;
	double a12 = 2 * (C1.y - B2.y);
	double b12 = 2 * (C1.x - B2.x);
	double c12 = C1.y * C1.y + C1.x * C1.x - B2.y * B2.y - B2.x * B2.x;

	//for A2B1C2
	double a21 = 2 * (B1.y - A2.y);
	double b21 = 2 * (B1.x - A2.x);
	double c21 = B1.y * B1.y + B1.x * B1.x - A2.y * A2.y - A2.x * A2.x;
	double a22 = 2 * (C2.y - B1.y);
	double b22 = 2 * (C2.x - B1.x);
	double c22 = C2.y * C2.y + C2.x * C2.x - B1.y * B1.y - B1.x * B1.x;
	//centers
	Point2f O1, O2;
	O1.x = (a11 * c12 - a12 * c11) / (a11 * b12 - a12 * b11);
	O1.y = (c11 * b12 - c12 * b11) / (a11 * b12 - a12 * b11);

	O2.x = (a21 * c22 - a22 * c21) / (a21 * b22 - a22 * b21);
	O2.y = (c21 * b22 - c22 * b21) / (a21 * b22 - a22 * b21);

	// S_A1C1A2
	double R_A1C1A2;
	Point2f O_A1C1A2;
	comCirCenterRadius(A1, C1, A2, &R_A1C1A2, &O_A1C1A2);

	// S_A1C1C2
	double R_A1C1C2;
	Point2f O_A1C1C2;
	comCirCenterRadius(A1, C1, C2, &R_A1C1C2, &O_A1C1C2);

	// S_A2C2A1
	double R_A2C2A1;
	Point2f O_A2C2A1;
	comCirCenterRadius(A2, C2, A1, &R_A2C2A1, &O_A2C2A1);

	// S_A2C2C1
	double R_A2C2C1;
	Point2f O_A2C2C1;
	comCirCenterRadius(A2, C2, C1, &R_A2C2C1, &O_A2C2C1);

	//using median computes centers and radius
	vector<double> tempR, tempOX, tempOY;
	// R 
	tempR.push_back(R1);
	tempR.push_back(R2);
	tempR.push_back(R_A1C1A2);
	tempR.push_back(R_A1C1C2);
	tempR.push_back(R_A2C2A1);
	tempR.push_back(R_A2C2C1);

	sort(tempR.begin(), tempR.end());
	*estimateR = (tempR[2] + tempR[3]) / 2.0;
	// O.x
	tempOX.push_back(O1.x);
	tempOX.push_back(O2.x);
	tempOX.push_back(O_A1C1A2.x);
	tempOX.push_back(O_A1C1C2.x);
	tempOX.push_back(O_A2C2A1.x);
	tempOX.push_back(O_A2C2C1.x);

	sort(tempOX.begin(), tempOX.end());
	double estimateOX = (tempOX[2] + tempOX[3]) / 2.0;
	//O.y
	tempOY.push_back(O1.y);
	tempOY.push_back(O2.y);
	tempOY.push_back(O_A1C1A2.y);
	tempOY.push_back(O_A1C1C2.y);
	tempOY.push_back(O_A2C2A1.y);
	tempOY.push_back(O_A2C2C1.y);

	sort(tempOY.begin(), tempOY.end());
	double estimateOY = (tempOY[2] + tempOY[3]) / 2.0;
	*estimateO = Point2f(estimateOX, estimateOY);
	return true;

}


Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er = 1.e-6)
{
	// SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin, Eigen::ComputeThinU | Eigen::ComputeThinV);
	// the result of SVD: UVD
	Eigen::MatrixXd U = svd_holder.matrixU();
	Eigen::MatrixXd V = svd_holder.matrixV();
	Eigen::MatrixXd D = svd_holder.singularValues();

	// construct S
	Eigen::MatrixXd S(V.cols(), U.cols());
	S.setZero();
	// inverse singularValues  and store in S
	for (unsigned int i = 0; i < D.size(); ++i)
	{

		if (D(i, 0) > er) {
			S(i, i) = 1 / D(i, 0);
		}
		else {
			S(i, i) = 0;
		}
	}

	// pinv_matrix = V * S * U^T
	return V * S * U.transpose();
}
Vec3d  refine_center_radius(vector<Point>& circle_pt, Vec3f& center_radius)
{
	int pt_num = circle_pt.size();
	MatrixXd A(pt_num, 3);
	MatrixXd E(pt_num, 1);
	Point2f circle_center(center_radius[0], center_radius[1]);
	double circle_r = center_radius[2];
	// construct coefficient matrix A and value matrix E
	for (int i = 0; i < pt_num; i++)
	{
		A(i, 0) = circle_pt[i].x - circle_center.x;
		A(i, 1) = circle_pt[i].y - circle_center.y;
		A(i, 2) = circle_r;
		E(i) = A(i, 0) * A(i, 0) + A(i, 1) * A(i, 1) - A(i, 2) * A(i, 2);
	}
	// pseudoinverse of A: use SVD and default error is 0
	Vec3d center_radius_refinement;
	center_radius_refinement[0] = center_radius[0] + (0.5 * pinv_eigen_based(A) * E)(0);
	center_radius_refinement[1] = center_radius[1] + (0.5 * pinv_eigen_based(A) * E)(1);
	center_radius_refinement[2] = center_radius[2] + (0.5 * pinv_eigen_based(A) * E)(2);
	return center_radius_refinement;
}

void CircleDetect::coCircleGroupArcs(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& out_groupedArcs,
	std::vector<std::vector<cv::Point>>& out_groupedArcsThreePt, std::vector<cv::Vec3f>& out_recordOR)
{
	vector<Point> vec;
	Point temp = Point(0, 0);
	vec.push_back(temp);
	for (int i = 0; i < in_edgelists.size(); i++)
	{
		int leng = in_edgelists[i].size();
		if (leng == 1) { continue; }//长度为1的段去除
		// put edge1 into outEdgeList first
		vector<Point> CirPt;
		vector<vector<Point> >outEdgeList;
		outEdgeList.push_back(in_edgelists[i]);

		Vec3f groupedOR;
		vector<Point> outThreePt;
		Point start = in_edgelists[i].front();
		Point end = in_edgelists[i].back();
		Point mid = in_edgelists[i][leng / 2];
		outThreePt.push_back(start);
		outThreePt.push_back(end);
		outThreePt.push_back(mid);//记录边缘段的起点、中点、终点

		vector<Vec3f> CenterRadius;//记录圆弧的半径和圆心
		//iterate to find the grouping arcs找到一个组的圆弧
		for (int j = 0; j < in_edgelists.size(); j++)
		{

			if (j == i || in_edgelists[j].size() == 1) { continue; }
			else
			{
				// group edgelist[j] with outEdgeList
				bool flag = false;
				int pass = 0;// record whether two edges can be grouped
				Vec3f temp1_center_radius;// record the two radii and centers for paired two arcs
				Vec3f temp2_center_radius;// temp1,temp2
				for (int k = 0; k < outEdgeList.size(); k++)
				{
					//两个弧是否属于同一个圆：是则返回圆的中心和半径，和两个弧的圆的中心和半径
					twoArcsCenterRadius(outEdgeList[k], in_edgelists[j], &flag, &temp1_center_radius, &temp2_center_radius, m_params.centerDis, m_params.radiusDis);
					if (flag)//can be grouped  
					{
						pass++;
						CenterRadius.push_back(temp1_center_radius);// record two radii and centers for paired arcs
						CenterRadius.push_back(temp2_center_radius);
					}
				}//endfor
				//
				if (pass == outEdgeList.size())// edgelist[j] can be grouped with  all before outEdgeList
				{
					outEdgeList.push_back(in_edgelists[j]);
					Point start2 = in_edgelists[j].front();
					Point end2 = in_edgelists[j].back();
					Point mid2 = in_edgelists[j][in_edgelists[j].size() / 2];
					outThreePt.push_back(start2);
					outThreePt.push_back(end2);
					outThreePt.push_back(mid2);
					in_edgelists[j] = vec;//归零
				}//endif
			}//endelse
		}//inner endfor
		// put points together 
		for (int l1 = 0; l1 < outEdgeList.size(); l1++)
		{
			for (int l2 = 0; l2 < outEdgeList[l1].size(); l2++)
			{
				CirPt.push_back(outEdgeList[l1][l2]);
			}
		}
		// estimate center and radius for single arc
		if (outEdgeList.size() == 1)// edgelist[i] didn't find grouped arcs
		{
			double singleR;
			Point2f singleO;
			estimateSingleCenterRadius(in_edgelists[i], &singleR, &singleO);
			if (std::isnan(singleO.y) || std::isnan(singleO.x) || std::isnan(singleR))
			{
				continue;
			}
			groupedOR = Vec3f(singleO.x, singleO.y, singleR);
		}
		// estimate center and radius for two arcs
		if (outEdgeList.size() == 2)// edgelist[i] finds one arc to pair
		{
			double TwoR;
			Point2f TwoO;
			estimateCenterRadius(outEdgeList[0], outEdgeList[1], &TwoR, &TwoO);
			if (std::isnan(TwoO.y) || std::isnan(TwoO.x) || std::isnan(TwoR))
			{
				continue;
			}
			groupedOR = Vec3f(TwoO.x, TwoO.y, TwoR);
		}
		//estimate center and radius for three or more arcs
		if (outEdgeList.size() > 2)// edgelist[i] finds more than two arcs, then take the median as estimated value
		{
			vector<double> temp_r;
			vector<float> temp_x;
			vector<float> temp_y;
			for (int i = 0; i < CenterRadius.size(); i++)
			{
				temp_x.push_back(CenterRadius[i][0]);
				temp_y.push_back(CenterRadius[i][1]);
				temp_r.push_back(CenterRadius[i][2]);
			}
			sort(temp_r.begin(), temp_r.end());
			sort(temp_x.begin(), temp_x.end());
			sort(temp_y.begin(), temp_y.end());
			// find the median for more than two arcs
			double MoreR;
			Point2f MoreO;
			int num_center_radius = CenterRadius.size();
			if (num_center_radius / 2 == 0)// even number
			{
				MoreR = (temp_r[num_center_radius / 2] + temp_r[(num_center_radius / 2) + 1]) / 2.0;
				MoreO.x = (temp_x[num_center_radius / 2] + temp_x[(num_center_radius / 2) + 1]) / 2;
				MoreO.y = (temp_y[num_center_radius / 2] + temp_y[(num_center_radius / 2) + 1]) / 2;

			}
			else {
				MoreR = temp_r[(num_center_radius + 1) / 2.0];
				MoreO.x = temp_x[(num_center_radius + 1) / 2.0];
				MoreO.y = temp_y[(num_center_radius + 1) / 2.0];

			}
			if (std::isnan(MoreO.y) || std::isnan(MoreO.x) || std::isnan(MoreR))
			{
				continue;
			}
			groupedOR = Vec3f(MoreO.x, MoreO.y, MoreR);
		}
		// refine the estimated center and radius优化圆心和半径

		Vec3d center_radius_refinement = refine_center_radius(CirPt, groupedOR);
		// finally detected grouped arcs, three points on it, and recoreded center and radius. 
		in_edgelists[i] = vec;
		out_groupedArcs.push_back(CirPt);
		out_groupedArcsThreePt.push_back(outThreePt);
		out_recordOR.push_back(center_radius_refinement);
	}//endfor outer
}


// 10 圆验证
bool circleVerify(vector<double>& x, vector<double>& y, int N, vector<Point>& stEdMid, Point2f O, double R, double* pe, double* angle)
{
	*pe = 0;
	*angle = 0;
	double pxc = O.x;
	double pyc = O.y;


	if (isinf(R) || isnan(R) || (R <= 0))
	{
		return false;
	}
	else
	{
		/*------------compute inliers after fitting------------*/
		/* first compute the distance from the edge point to the fitted circle, the distance constraint is
		then compute the angle between point normals and circle normals,
		if these two constraints hold, further compute the inlier ratio and the span angle.
		*/

		//first compute circular points' normals using the contour method
		vector<Point2f> DT;
		vector<Point2f> D1(N), D2(N);
		for (int j = 0; j < N - 1; j++)
		{   // differences between two consecutive points
			float L = sqrt((x[j] - x[j + 1]) * (x[j] - x[j + 1]) + (y[j] - y[j + 1]) * (y[j] - y[j + 1]));
			D1[j].x = (x[j] - x[j + 1]) / L;
			D1[j].y = (y[j] - y[j + 1]) / L;
			D2[j + 1].x = (x[j] - x[j + 1]) / L;
			D2[j + 1].y = (y[j] - y[j + 1]) / L;
		}//endfor
		float LL = sqrt((x.back() - x.front()) * (x.back() - x.front()) + (y.back() - y.front()) * (y.back() - y.front()));
		Point2f end2start;
		end2start.x = (x.back() - x.front()) / LL;
		end2start.y = (y.back() - y.front()) / LL;
		D1.push_back(end2start);
		D2.front() = end2start;
		//D=D1+D2 Normal
		vector<Point2f> D(N), Nom(N);
		for (int i = 0; i < N; i++)
		{
			D[i].x = D1[i].x + D2[i].x;
			D[i].y = D1[i].y + D2[i].y;
			double LL1 = sqrt(D[i].x * D[i].x + D[i].y * D[i].y);
			Nom[i].x = -D[i].y / LL1;
			Nom[i].y = D[i].x / LL1;
		}//endfor


		//compute circular points' normals using fitted centers(equation)
		int inlierNum = 0;
		for (int i = 0; i < N; i++) {
			double dx = x[i] - pxc;

			double dy = y[i] - pyc;
			//point normals connecting the center
			double distP2O = sqrt(dx * dx + dy * dy);
			double d = abs(distP2O - R);
			if (d <= 2)// distance less than 2 pixel
			{

				float theta = atan2(dy, dx);
				Point2f circleNormal = Point2f(cos(theta), sin(theta));

				float cosPtCirNormal2 = (Nom[i].x * circleNormal.x + Nom[i].y * circleNormal.y);
				if (abs(cosPtCirNormal2) >= cos(22.5 / 180 * CV_PI))// angle<=22.5°remember float calculations
				{
					inlierNum++;
				}//endif

			}//endif
			//waitKey();
		} //end-for

		double inlierEdgeRatio = (double)inlierNum / N;
		double inlierRatio = 0, spanAngle = 0;

		if (inlierEdgeRatio >= 0.2)// take up 0.5\0.8 of edge points
		{

			inlierRatio = inlierNum / (2 * CV_PI * R);

		}//endif
		*pe = inlierRatio;
		*angle = spanAngle;
	}//endelse
	return true;
}

vector<CircleData> circleEstimateGroupedArcs(vector<vector<Point>>& groupedArcs, vector<Vec3f>& recordOR, vector<vector<Point>>& groupedArcsThreePt, float& T_inlier)
{
	vector<CircleData> addCircles;
	for (int i = 0; i < groupedArcs.size(); i++)
	{
		CircleData fitCircle;
		vector<double> X, Y;// point coordinates
		vector<Point> stEdMid;// arcs start, mid and end points;
		stEdMid = groupedArcsThreePt[i];
		double groupedR = recordOR[i][2];
		Point2f groupedO(recordOR[i][0], recordOR[i][1]);
		for (int j = 0; j < groupedArcs[i].size(); j++)
		{	/* or for (auto j = groupedArcs[i].begin(); j != groupedArcs[i].end(); j++)*/
			Y.push_back(groupedArcs[i][j].y);// or X.push_back((*j).x)
			X.push_back(groupedArcs[i][j].x);

		}//endfor

		//fit
		double  inlierRatio, spanAngle;
		circleVerify(X, Y, X.size(), stEdMid, groupedO, groupedR, &inlierRatio, &spanAngle);

		// inlier verification
		if (inlierRatio >= T_inlier && groupedR >= 5)//&& spanAngle>=1 && spanAngle<=0 && spanAngle >= T_angle spanned angle can be omitted
		{
			fitCircle.center = Point2f(groupedO.x, groupedO.y);
			fitCircle.radius = groupedR;
			fitCircle.ratio = inlierRatio;
			addCircles.push_back(fitCircle);

		}//endif

	}//endfor                                                                      
	return addCircles;
}

bool estimateClosedCenterRadius(vector<Point>& A1B1C1, double* estimateR, Point2f* estimateO)
{
	// sample 5pts evenly
	int ptNum = A1B1C1.size();
	Point A = A1B1C1[ptNum / 5];
	Point B = A1B1C1[2 * ptNum / 5];
	Point C = A1B1C1[3 * ptNum / 5];
	Point D = A1B1C1[4 * ptNum / 5];
	Point E = A1B1C1[ptNum - 7];//A1B1C1[ptNum-10]



	// five new inscribed triangles
	//ACD
	double R_ACD;
	Point2f O_ACD;
	comCirCenterRadius(A, C, D, &R_ACD, &O_ACD);

	//ACE
	double R_ACE;
	Point2f O_ACE;
	comCirCenterRadius(A, C, E, &R_ACE, &O_ACE);


	//ABD
	double R_ABD;
	Point2f O_ABD;
	comCirCenterRadius(A, B, D, &R_ABD, &O_ABD);


	//BCE
	double R_BCE;
	Point2f O_BCE;
	comCirCenterRadius(B, C, E, &R_BCE, &O_BCE);


	//DBE
	double R_DBE;
	Point2f O_DBE;
	comCirCenterRadius(D, B, E, &R_DBE, &O_DBE);


	// using median computes centers and radius
	vector<double> tempR, tempOX, tempOY;
	// R 
	tempR.push_back(R_ACD);
	tempR.push_back(R_ACE);
	tempR.push_back(R_ABD);
	tempR.push_back(R_BCE);
	tempR.push_back(R_DBE);
	//OX
	tempOX.push_back(O_ACD.x);
	tempOX.push_back(O_ACE.x);
	tempOX.push_back(O_ABD.x);
	tempOX.push_back(O_BCE.x);
	tempOX.push_back(O_DBE.x);

	//OY
	tempOY.push_back(O_ACD.y);
	tempOY.push_back(O_ACE.y);
	tempOY.push_back(O_ABD.y);
	tempOY.push_back(O_BCE.y);
	tempOY.push_back(O_DBE.y);


	sort(tempR.begin(), tempR.end());
	*estimateR = tempR[2];

	sort(tempOX.begin(), tempOX.end());
	double estimateOX = tempOX[2];

	sort(tempOY.begin(), tempOY.end());
	double estimateOY = tempOY[2];
	*estimateO = Point2f(estimateOX, estimateOY);

	return true;

}
vector<CircleData> circleEstimateClosedArcs(vector<vector<Point>>& closedArcs, float& T_inlier_closed)
{


	Mat pre = Mat(500, 700, CV_8UC3, Scalar(255, 255, 255));


	vector<CircleData> addCircles;

	for (int i = 0; i < closedArcs.size(); i++)
	{

		CircleData fitCircle;
		vector<double> X, Y;
		vector<Point> threePt;
		double closedR;
		Point2f closedO;
		estimateClosedCenterRadius(closedArcs[i], &closedR, &closedO);// estimate the center and radius
		//cout << i << "ClosedCenter: " << closedO.x << " " << closedO.y << endl;

		int r = rand() % 256;
		int g = rand() % 256;
		int b = rand() % 256;
		Scalar colorClosedEdgesPre = Scalar(b, g, r);
		for (int j = 0; j < closedArcs[i].size(); j++)
		{	/* or for (auto j = groupedArcs[i].begin(); j != groupedArcs[i].end(); j++)*/
			Y.push_back(closedArcs[i][j].y);// or X.push_back((*j).x)
			X.push_back(closedArcs[i][j].x);
			circle(pre, closedArcs[i][j], 1, colorClosedEdgesPre);
		}//endfor

		//fit
		double inlierRatio, spanAngle;
		circleVerify(X, Y, X.size(), threePt, closedO, closedR, &inlierRatio, &spanAngle);


		// inlier verification

		if (inlierRatio >= T_inlier_closed && closedR >= 4)//inlierRatio >=0.5 
		{
			fitCircle.center = Point2f(closedO.x, closedO.y);
			fitCircle.radius = closedR;
			fitCircle.ratio = inlierRatio;
			addCircles.push_back(fitCircle);

		}//endif
		//cout << "Add circle done" << endl;
	}//endfor
	return addCircles;
}

void CircleDetect::toalcircle(std::vector<std::vector<cv::Point>>& in_edgelists, std::vector<std::vector<cv::Point>>& in_groupedArcs,
	std::vector<std::vector<cv::Point>>& in_groupedArcsThreePt,
	std::vector<cv::Vec3f>& in_recordOR, std::vector<CircleData>& out_Circles)
{
	vector<CircleData> groupedCircles;// grouped arcs
	groupedCircles = circleEstimateGroupedArcs(in_groupedArcs, in_recordOR, in_groupedArcsThreePt, m_params.inlier);//fit grouped arcs

	// closed arcs闭合链
	vector<CircleData> closedCircles;// closedCircles
	//闭合圆拟合
	closedCircles = circleEstimateClosedArcs(in_edgelists, m_params.closedInlier);// fit closed edges
	//put grouped and closed circles together
	if (!groupedCircles.empty())
	{
		out_Circles = groupedCircles;
	}
	if (!closedCircles.empty())
	{
		for (auto it = closedCircles.begin(); it != closedCircles.end(); it++)
		{
			out_Circles.push_back(*it);
		}
	}
}



// 11 去除相似圆
void CircleDetect::clusterCircles(std::vector<CircleData>& in_Circles, std::vector<CircleData>& out_Circles)
{
	while (!in_Circles.empty())
	{
		vector<CircleData> simCircles;// similar circles
		CircleData circle1 = in_Circles.front();
		for (auto it = in_Circles.begin(); it != in_Circles.end();)
		{
			CircleData circle2 = *it;
			//判断圆是否相似
			double disCircles = sqrt((circle1.center.x - circle2.center.x) * (circle1.center.x - circle2.center.x) +
				(circle1.center.y - circle2.center.y) * (circle1.center.y - circle2.center.y) + 
				(circle1.radius - circle2.radius) * (circle1.radius - circle2.radius));
			if ((disCircles <= 5))// cluster 5\10
			{
				simCircles.push_back(circle2);// put together
				it = in_Circles.erase(it);// delete it
			}
			else { it++; }
		}
		// find the maximum inlier ratio in simCircles
		//在相似圆中找到圆率最大的
		if (simCircles.size() > 1)
		{
			sort(simCircles.begin(), simCircles.end(), [](const CircleData& a, const CircleData& b){return a.ratio > b.ratio;});
		}// sort default: increase order; cmpInlier is decrease order
		out_Circles.push_back(simCircles.front());
	}
}

void showPoints(cv::Mat gray, std::vector<std::vector<cv::Point>>& in_edgelists)
{
	namedWindow("img", 1);
	Mat show;
	cvtColor(gray, show, COLOR_GRAY2RGB);
	for (auto& points : in_edgelists)
	{
		for (auto& point : points)
		{
			int x = point.x;
			int y = point.y;
			//circle(show, Point(x,y), 1, { 0, 255, 0 }, 1);			
			show.at<Vec3b>(y, x) = { 0, 255, 0 };
		}

	}
	imshow("img", show);
	imwrite("img.jpg", show);
}

void CircleDetect::detect(cv::Mat image)
{
	// 0 初始化
	if (image.empty())
	{
		std::cout << "image is empty!" << endl;
		return;
	}
	switch (image.channels())
	{
	case 1:
		this->m_gray = image;
		break;
	case 3:
		cvtColor(image, this->m_gray, cv::COLOR_BGR2GRAY);
		break;
	case 4:
		cvtColor(image, this->m_gray, cv::COLOR_BGRA2GRAY);
		break;
	default:
		break;
	}


	// 1 滤波 并 检测边缘段
	GaussianBlur(m_gray, m_gray, Size(m_params.filtSize, m_params.filtSize), 0, 0);
	EDPF testEDPF = EDPF(m_gray);
	std::vector<std::vector<cv::Point>> t_segments = testEDPF.getSegments();
	vector<vector<Point>> t_edgeList;
	for (int i = 0; i < t_segments.size(); i++)
	{
		if (t_segments[i].size() >= m_params.minLen)
		{
			t_edgeList.push_back(t_segments[i]);
		}
	}
	if (t_segments.size() == 0)
	{
		m_results = Mat(1, 1, CV_32FC1, Scalar(0));
		return;
	}

	// 2 检测闭合边缘段
	vector<vector<Point>> t_closedSegments;
	extractClosedEdges(t_edgeList, t_closedSegments);

	// 3 以直代曲
	vector<vector<Point>> t_rdpSegments;
	RDP(t_edgeList, t_rdpSegments);

	// 4 链中角度较大的拐点处断开
	vector<vector<Point>> t_saSegments;
	vector<vector<Point>> t_saEdgelists;
	rejectSharpTurn(t_rdpSegments, t_edgeList, t_saSegments, t_saEdgelists);


	// 5 直线弧线交界处断开
	vector<vector<Point>> t_dlSegments;
	vector<vector<Point>> t_dlsaEdgelists;
	detectline(t_saSegments, t_saEdgelists, t_dlSegments, t_dlsaEdgelists);


	// 6 根据链的走向、弯曲方向变化断开链
	vector<vector<Point>> t_dfSegments;
	vector<vector<Point>> t_dfsaEdgelists;
	detectInflexPt(t_dlSegments, t_dlsaEdgelists, t_dfSegments, t_dfsaEdgelists);

	// 7 删除段的链和直线
	detectshort(t_dfsaEdgelists);
	
	// 8 提取闭合边和非闭合边
	vector<vector<Point>> t_cSegments;
	vector<vector<Point>> t_ncSegments;
	extractClosedEdges2(t_dfsaEdgelists, t_ncSegments, t_cSegments);

	// 9 非闭合曲线分组
	vector<vector<Point> > groupedArcs;
	vector<vector<Point> > groupedArcsThreePt;
	vector<Vec3f>recordOR;
	vector<CircleData> groupedCircles;
	coCircleGroupArcs(t_ncSegments, groupedArcs, groupedArcsThreePt, recordOR);

	// 10 圆验证
	std::vector<CircleData> totalCircles;
	toalcircle(t_cSegments, groupedArcs, groupedArcsThreePt, recordOR, totalCircles);

	// 11 去除相似圆
	std::vector<CircleData> resultCircles;
	clusterCircles(totalCircles, resultCircles);

	if (!resultCircles.empty())
	{
		Mat result1((int)resultCircles.size(), 4, CV_32FC1, Scalar(0));
		for (int i = 0; i < (int)resultCircles.size(); i++)
		{
			float* img_ptr = (float*)result1.ptr(i);
			img_ptr[0] = resultCircles[i].center.x;
			img_ptr[1] = resultCircles[i].center.y;
			img_ptr[2] = resultCircles[i].radius;
			img_ptr[3] = resultCircles[i].ratio;
		}
		m_results = result1.clone();
	}
	else
	{
		Mat result1(1, 1, CV_32FC1, Scalar(0));
		m_results = result1.clone();
	}

	//showPoints(m_gray, t_ncSegments);
}