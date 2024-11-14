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
  
![6](https://github.com/user-attachments/assets/c8f18057-257d-42c2-b2db-8d548e775a99)
