from DVisionCircleDetect import DVisionCircleDetect

# 创建模型
class CircleDetect():
    def __new__(self):
        return DVisionCircleDetect.CircleDetect()

    def setSplit(self, lineLen, sharpAngle, SegmentCurvature):
        """
            .   @brief 输入边缘段分割的参数
            .   @param lineLen: 直线弧线交界处分割，大于lineLen的认为是直线 【160】
            .   @param sharpAngle: 线段拐点处分割，sharpAngle 角度大于阈值处分割 【90】
            .   @param SegmentCurvature: 弧线段曲率，根据弯曲程度断开 【0.01】
            """
        pass

    def setInlier(self, inlier, closedInlier):
        """
            .   @brief 输入内点率
            .   @param inlier: 非闭合圆的内点率
            .   @param closedInlier: 闭合圆的内点率
            """
        pass

    def setDistance(self, centerDis, radiusDis):
        """
            .   @brief 输入判断相似圆的距离参数
            .   @param centerDis: 圆中心的距离
            .   @param radiusDis: 圆半径之差
            """
        pass

    def setFiltSize(self, filtSize):
        """
            .   @brief 滤波尺寸
            .   @param filtSize: 滤波尺寸
            """
        pass

    def setMinLen(self, minLen):
        """
            .   @brief 输入最小线段长度
            .   @param minLen: 是最小线段长度
            """
        pass

    def detect(self, image):
        """
            .   @brief 检测
            .   @param image: 图像矩阵
            """
        pass

    def getResult(self):
        """
            .   @brief 返回结果
            .   @param out: 图像矩阵,大小为 N*4 ，N为圆的数量，0~3：圆心x,y，半径，内点率
            """
        pass



