import cv2
import numpy as np
import random
import DVisionCircleDetect
import math




# 获取图像
imageID = 2
if imageID == 1:
    src = cv2.imread('image/01.jpg')
if imageID == 2:
    src = cv2.imread('image/2.jpg')


# 创建模型
# 创建模型类
m_CircleDetect = DVisionCircleDetect.CircleDetect()
# 设置参数
m_CircleDetect.setFiltSize(5)
m_CircleDetect.setMinLen(32)
m_CircleDetect.setDistance(5, 5)
m_CircleDetect.setInlier(0.5, 0.5)
m_CircleDetect.setSplit(160, 60, 0.01)
# 检测
m_CircleDetect.detect(src)
# 返回结果
result = m_CircleDetect.getResult()

(dstrows, dstcols) = (result.shape[0], result.shape[1])
print(dstrows, dstcols)

# 显示结果图像
radius = 3
color = (0, 0, 255)
thickness = 1  # 填充实心圆
radius2 = 1
color2 = (0, 255, 0)
thickness2 = 2  # 填充实心圆
if dstcols > 1:
    for i in range(0, dstrows):
        cv2.circle(src, (round(result[i][0]), round(result[i][1])), round(result[i][2]), color, thickness)
        cv2.circle(src, (round(result[i][0]), round(result[i][1])), radius2, color2, thickness2)
        cv2.putText(src, "ratio: " + str("{:.3f}".format(result[i][3])), (round(result[i][0])+5, round(result[i][1])), cv2.FONT_HERSHEY_SIMPLEX,
                0.5, (0, 0, 255), 1, cv2.LINE_AA)
    cv2.imshow("alllines", src)
    cv2.imwrite('123.jpg',src)
    cv2.waitKey()

