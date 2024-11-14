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

3.shape match轮廓匹配<br>
![轮廓匹配](图像算法库说明/img/9.jpg)<br>
![轮廓匹配](图像算法库说明/img/10.jpg)<br>
![轮廓匹配](图像算法库说明/img/12.jpg)<br>


4.形状/灰度匹配<br>
![形状匹配](图像算法库说明/img/19.jpg)<br>

5.棋盘格标定<br>
![棋盘格标定](图像算法库说明/img/20.jpg)<br>

5.形状检测-圆<br>
![形状检测](图像算法库说明/img/17.jpg)<br>


6.边缘检测<br>
![边缘检测](图像算法库说明/img/18.jpg)<br>

# 其他参考
https://github.com/DennisLiu1993/Fastest_Image_Pattern_Matching
https://github.com/meiqua/shape_based_matching
https://github.com/sdg002/RANSAC
https://github.com/glassechidna/zxing-cpp
