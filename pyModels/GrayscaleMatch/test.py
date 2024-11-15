import cv2
import numpy as np
import random
import DVisionGrayscaleMatch
import math


# 获取图像
imageID = 0
if imageID == 0:
    dst = cv2.imread('image/2temp_image/img1.bmp')
    src = cv2.imread('image/1src_image/img1.bmp')
if imageID == 1:
    dst = cv2.imread('image/6_12.png')
    src = cv2.imread('image/6_13.bmp')
if imageID == 2:
    dst = cv2.imread('image/641.jpg')
    src = cv2.imread('image/640.png')
if imageID == 3:
    dst = cv2.imread('image/b2.png')
    src = cv2.imread('image/a1.bmp')
if imageID == 4:
    dst = cv2.imread('image/test11.png')
    src = cv2.imread('image/8du.jpg')
if imageID == 5:
    dst = cv2.imread('image/tep1.bmp')
    src = cv2.imread('image/2D00.png')

# 创建模型
# 创建模型类
m_GrayscaleMatch = DVisionGrayscaleMatch.GrayscaleMatch()
# 设置参数
m_GrayscaleMatch.setTempimage(dst)
m_GrayscaleMatch.setScore(0.5)
m_GrayscaleMatch.setAngle(0,360)
# m_GrayscaleMatch.setPyramidLayer(16)
# 匹配
m_GrayscaleMatch.match(src)
# 返回结果
result = m_GrayscaleMatch.getResult()


(dstrows, dstcols) = (result.shape[0], result.shape[1])
print(dstrows, dstcols)

# 显示结果图像
radius = 3
color = (0, 0, 255)
thickness = 3  # 填充实心圆
radius2 = 1
color2 = (0, 255, 0)
thickness2 = 1  # 填充实心圆
if dstcols > 1:
    for i in range(0, dstrows):
        cv2.circle(src, (round(result[i][0]), round(result[i][1])), radius, color, thickness)

        cv2.putText(src, "s: " + str("{:.3f}".format(result[i][2])), (round(result[i][0])+5, round(result[i][1])), cv2.FONT_HERSHEY_SIMPLEX,
                0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "a: " + str("{:.3f}".format(result[i][3])), (round(result[i][0] + 5), round(result[i][1]+16)), cv2.FONT_HERSHEY_SIMPLEX, 0.5,
                (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "x: " + str("{:.3f}".format(result[i][0])), (round(result[i][0] + 5), round(result[i][1] + 32)), cv2.FONT_HERSHEY_SIMPLEX,
                0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "y: " + str("{:.3f}".format(result[i][1])), (round(result[i][0]) + 5, round(result[i][1] + 48)), cv2.FONT_HERSHEY_SIMPLEX,
                0.5, (0, 0, 255), 1, cv2.LINE_AA)
    cv2.imshow("alllines", src)
    cv2.imwrite("a.jpg", src)
    cv2.waitKey()

