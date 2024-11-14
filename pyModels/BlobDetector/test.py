import cv2
from DVisionBlobDetector import DVisionBlobDetector
import numpy as np
# 获取图像
imageID = 2
if imageID == 1:
    src = cv2.imread('image/640.png')
if imageID == 2:
    src = cv2.imread('image/a1.bmp')
if imageID == 3:
    src = cv2.imread('image/2DMeasure.bmp')

# 创建类
m_BlobDetector = DVisionBlobDetector.BlobDetector()
# 设置参数
#  排序
m_BlobDetector.setSort(1)
# 阈值
m_BlobDetector.setThreshold(1, 70, 255)
#  连通域颜色
m_BlobDetector.setColor(1)
#  面积筛选
m_BlobDetector.setArea(1, 2500)
#  周长筛选
m_BlobDetector.setPerimeter(0, 50)
#  轴比筛选
m_BlobDetector.setInertiaRatio(0, 0.1, 1)
#  矩形度筛选
m_BlobDetector.setRectangularity(0, 0.9, 1)
# 圆形度筛选
m_BlobDetector.setCircularity(0, 0.8, 1)
#  质心距离筛选
m_BlobDetector.setCentroidOffset(0, 2)
#  最大数量
m_BlobDetector.setmaxTargets(0)


# 检测
m_BlobDetector.detect(src)

# 获取匹配结果
result = m_BlobDetector.getResults()
contours = m_BlobDetector.getContours()
externalContours = m_BlobDetector.getExternalContours()
minRects = m_BlobDetector.getMinRects()
binary = m_BlobDetector.getBinaryimage()
rbinary = m_BlobDetector.getResultBinary()
(dstrows, dstcols) = (result.shape[0], result.shape[1])
print(dstrows)

cv2.imshow("binary", binary)
cv2.imshow("rbinary", rbinary)

radius1 = 1
color1 = (0, 0, 255)
thickness = 2  # 填充实心圆
radius2 = 1
color2 = (0, 255, 0)
thickness2 = 1  # 填充实心圆
if dstcols > 1:
    for i in range(0, dstrows):
        # 最小外接矩形
        cv2.line(src,(round(minRects[i][5]), round(minRects[i][6])), (round(minRects[i][7]), round(minRects[i][8])),
                 color2, thickness2)
        cv2.line(src, (round(minRects[i][7]), round(minRects[i][8])), (round(minRects[i][9]), round(minRects[i][10])),
                 color2, thickness2)
        cv2.line(src, (round(minRects[i][9]), round(minRects[i][10])), (round(minRects[i][11]), round(minRects[i][12])),
                 color2, thickness2)
        cv2.line(src, (round(minRects[i][11]), round(minRects[i][12])), (round(minRects[i][5]), round(minRects[i][6])),
                 color2, thickness2)

        # 轮廓点
        for j in range(0, round(contours[i][0])):
            cv2.circle(src, (round(contours[i][2 * j + 1]), round(contours[i][2 * j + 2])), radius2, color2, thickness2)

        # 连通域参数
        cv2.circle(src, (round(result[i][0]), round(result[i][1])), radius1, color1, thickness)
        cv2.putText(src, "area: " + str("{:.3f}".format(result[i][2])), (round(result[i][0]) + 5, round(result[i][1])),
                    cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "perimeter: " + str("{:.3f}".format(result[i][3])),
                    (round(result[i][0]) + 5, round(result[i][1]) + 16), cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "circularity: " + str("{:.3f}".format(result[i][4])),
                    (round(result[i][0]) + 5, round(result[i][1]) + 32), cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "rectangularity: " + str("{:.3f}".format(result[i][5])),
                    (round(result[i][0]) + 5, round(result[i][1]) + 48), cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(src, "axleratio: " + str("{:.3f}".format(result[i][6])),
                    (round(result[i][0]) + 5, round(result[i][1]) + 64), cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, (0, 0, 255), 1, cv2.LINE_AA)

    cv2.imshow("src", src)
    cv2.waitKey()

