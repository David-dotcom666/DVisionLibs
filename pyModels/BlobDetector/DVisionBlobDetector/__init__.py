from DVisionBlobDetector import DVisionBlobDetector

# 创建模型
class BlobDetector():
    def __new__(self):
        return DVisionBlobDetector.BlobDetector()

    def setColor(self, blobColor):
        """
            .   @brief 设置连通域颜色
            .   @param blobColor: [连通域颜色] [类型：bool] [默认值：1] [0黑/1白]
            """
        pass
    def setmaxTargets(self, maxTargets):
        """
            .   @brief 设置连通域寻找最大数量
            .   @param maxTargets: [最大数量] [类型：int] [默认值：0] [0时找到所有连通域]
            """
        pass
    def setSort(self, sortMethod, sortFeature):
        """
            .   @brief 设置排序方式
            .   @param sortMethod: [排序方式] [类型：int] [默认值：1] [0：不排序，1：降序，2：升序]
            .   @param sortFeature: [排序特征] [类型：int] [默认值：0] [0：面积，1：周长，2：圆形度，3：矩形度，4：轴比，5：外接矩形中心x，6：中心y]
            """
        pass
    def setThreshold(self, autoThreshold, minThreshold, maxThreshold):
        """
            .   @brief 设置二值化阈值参数
            .   @param autoThreshold: [是否自动确定阈值] [类型：int] [默认值：1] [0：自定义阈值，1：双峰类型直方图自动确定阈值，2：单峰类型直方图动确定阈值]
            .   @param minThreshold: [低阈值] [类型：int] [默认值：155]
            .   @param maxThreshold: [高阈值] [类型：int] [默认值：255]
            """
        pass
    def setArea(self, filterByArea, minArea, maxArea):
        """
            .   @brief 设置面积筛选参数
            .   @param filterByArea: [是否使用面积筛选] [类型：bool] [默认值：1] [0：不使用，1：使用]
            .   @param minArea: [低阈值] [类型：float] [默认值：25]
            .   @param maxArea: [高阈值] [类型：float] [默认值：浮点数最大]
            """
        pass
    def setCircularity(self, filterByCircularity, minCircularity, maxCircularity):
        """
            .   @brief 设置圆形度筛选参数
            .   @param filterByCircularity: [是否使用圆形度筛选] [类型：bool] [默认值：0] [0：不使用，1：使用]
            .   @param minCircularity: [低阈值] [类型：float] [默认值：0.1]
            .   @param maxCircularity: [高阈值] [类型：float] [默认值：1]
            """
        pass
    def setInertiaRatio(self, filterByInertia, minInertiaRatio, maxInertiaRatio):
        """
            .   @brief 设置长短轴比筛选参数
            .   @param filterByCircularity: [是否使用长短轴比筛选] [类型：bool] [默认值：0] [0：不使用，1：使用]
            .   @param minCircularity: [低阈值] [类型：float] [默认值：0.1]
            .   @param maxCircularity: [高阈值] [类型：float] [默认值：1]
            """
        pass
    def setRectangularity(self, filterByRectangularity, minRectangularity, maxRectangularity):
        """
            .   @brief 设置长短轴比筛选参数
            .   @param filterByRectangularity: [是否使用长短轴比筛选] [类型：bool] [默认值：0] [0：不使用，1：使用]
            .   @param minRectangularity: [低阈值] [类型：float] [默认值：0.1]
            .   @param maxRectangularity: [高阈值] [类型：float] [默认值：1]
            """
        pass
    def setPerimeter(self, filterByPerimeter, minPerimeter, maxPerimeter):
        """
            .   @brief 设置周长筛选参数
            .   @param filterByArea: [是否使用周长筛选] [类型：bool] [默认值：0] [0：不使用，1：使用]
            .   @param minArea: [低阈值] [类型：float] [默认值：25]
            .   @param maxArea: [高阈值] [类型：float] [默认值：浮点数最大]
            """
        pass
    def setCentroidOffset(self, filterByCentroidOffset, minCentroidOffset, maxCentroidOffset):
        """
            .   @brief 设置质心偏移距筛选参数
            .   @param filterByArea: [是否使用质心偏移距筛选] [类型：bool] [默认值：0] [0：不使用，1：使用]
            .   @param minArea: [低阈值] [类型：float] [默认值：25]
            .   @param maxArea: [高阈值] [类型：float] [默认值：浮点数最大]
            """
        pass
    def detect(self, _image, _mask):
        """
            .   @brief 连通域检测函数
            .   @param _image: [检测图像] [类型：Mat]
            .   @param _mask: [掩码图像] [类型：Mat] [默认值：空]
            """
        pass
    def getResults(self):
        """
            .   @brief 获取检测结果，连通域一些参数
            .   @out: Mat类型
                      没有结果时为 1*1 大小值为0的矩阵；
                      有结果时为 N*8 大小的浮点数矩阵，N为连通域数量，每行参数为：[0,1: 质心的坐标x,y][2: 面积][3: 周长][4: 圆形度][5: 矩形度][6: 轴比][7: 质心偏移距]
            """
        pass
    def getContours(self):
        """
            .   @brief 获取检测结果，连通域所有轮廓点
            .   @out: Mat类型
                      没有结果时为 1*1 大小值为0的矩阵；
                      有结果时为 N*M 大小的浮点数矩阵，N为连通域数量，每行参数为：[0: 当前行对应连通域轮廓点的数量][1,2: 轮廓点的坐标x,y][3,4....]
            """
        pass
    def getExternalContours(self):
        """
            .   @brief 获取检测结果，连通域外轮廓点
            .   @out: Mat类型
                      没有结果时为 1*1 大小值为0的矩阵；
                      有结果时为 N*M 大小的浮点数矩阵，N为连通域数量，每行参数为：[0: 当前行对应连通域外轮廓点的数量][1,2: 轮廓点的坐标x,y][3,4....]
            """
        pass
    def getMinRects(self):
        """
            .   @brief 获取检测结果，连通域最小外接矩形
            .   @out: Mat类型
                      没有结果时为 1*1 大小值为0的矩阵；
                      有结果时为 N*13 大小的浮点数矩阵，N为连通域数量，每行参数为：[0,1: 中心点坐标][2: 角度][3,4: 短轴，长轴][5~12: 四个角点坐标]
            """
        pass
    def getBinaryimage(self):
        """
            .   @brief 获取检测结果，二值图像
            .   @out: Mat类型
            """
        pass
    def getResultBinary(self):
        """
            .   @brief 获取检测结果，连通域结果二值图像
            .   @out: Mat类型
            """
        pass