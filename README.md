# DVisionLibs
Traditional machine vision algorithm library, including code reading, blob detection, circle detection, line detection, shape matching, contour matching, chessboard calibration, etc

1.Read Barcode<br>
基于ZXing库,主要提供了非本地图像路径而是从变量进行识别的接口，支持C++和pyd接口<br>
用法：<br>
a.C++<br>
  首先按照说明配置路径,将项目属性的配置类型改成"应用程序.exe"<br>
  执行test.cpp<br>
b.python<br>
  执行pymodels下的test.py<br>
c.修改编译pyd<br>
  修改后注意将项目属性的配置类型改成"动态库.dll",将配置属性-高级-高级属性-目标文件扩展名改成:.pyd<br>
d.修改编译dll<br>
  修改后注意将项目属性的配置类型改成"动态库.dll",将配置属性-高级-高级属性-目标文件扩展名改成:.dll<br>
e.需要选择编译指定模块,在设置中设置生成库文件的名称,其他的模块需要在vs中设置选择从项目中排除
  
![barcode result](图像算法库说明/img/6.png)<br>

2.Blob Detection<br>
用法同上<br>
a.python版本中可以查看检测结果<br>
![blob result](图像算法库说明/img/7.png)<br>
![blob result](图像算法库说明/img/8.png)<br>

3.ShapeMatch/轮廓匹配<br>
![轮廓匹配](图像算法库说明/img/9.jpg)<br>
![轮廓匹配](图像算法库说明/img/10.jpg)<br>
![轮廓匹配](图像算法库说明/img/12.jpg)<br>


4.NCCMatch/GrayMatch灰度匹配<br>
	除了常用的分数,极性,角度范围,数量,重叠率,亚像素等还包括了支持mask功能,均值和标准差加速功能

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

![形状匹配](图像算法库说明/img/19.jpg)<br>

5.CameraCalibration/ChessboardCalibration棋盘格标定<br>
![棋盘格标定](图像算法库说明/img/20.jpg)<br>

6.CircleDetection/形状检测-圆<br>
![形状检测](图像算法库说明/img/17.jpg)<br>


7.EdgeDetection/边缘检测<br>
![边缘检测](图像算法库说明/img/18.jpg)<br>

# 其他参考
https://github.com/DennisLiu1993/Fastest_Image_Pattern_Matching<br>
https://github.com/meiqua/shape_based_matching<br>
https://github.com/sdg002/RANSAC<br>
https://github.com/glassechidna/zxing-cpp<br>
