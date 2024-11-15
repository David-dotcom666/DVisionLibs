from DVisionGrayscaleMatch import DVisionGrayscaleMatch

# 创建模型
class GrayscaleMatch():
    def __new__(self):
        return DVisionGrayscaleMatch.GrayscaleMatch()

    def setTempimage(self, temp, useMask, mask):
        """
            .   @brief 输入模板图像 和 掩码
            .   @param temp: 模板图像
            .   @param useMask: 是否使用mask 【bool】 【默认值：0】
            .   @param mask: 模板图像 【默认值：空】
            """
        pass

    def setScore(self, score):
        """
            .   @brief 输入匹配分数
            .   @param score: 匹配分数 【float】
            """
        pass

    def setMaxtargs(self, maxtargs):
        """
            .   @brief 输入最大匹配数量
            .   @param maxtargs: 最大匹配数量
            """
        pass

    def setPyramidLayer(self, pyramidLayer):
        """
            .   @brief 输入金字塔参数
            .   @param pyramidLayer: 金字塔参数
            """
        pass

    def setAngle(self, angleStart, angleRange):
        """
            .   @brief 输入匹配角度
            .   @param angleStart: 起始角度
            .   @param angleRange: 角度范围
            """
        pass

    def setMaxOverlap(self, setMaxOverlap):
        """
            .   @brief 输入最大重叠率
            .   @param setMaxOverlap: 最大重叠率
            """
        pass

    def setMean(self, useMean, mean):
        """
            .   @brief 使用亮度均值加速运算
            """
        pass

    def setSDV(self, useSDV, SDV):
        """
            .   @brief 使用对比度加速运算
            """
        pass

    def match(self, image):
        """
            .   @brief 匹配
            .   @param image: 匹配图像
            """
        pass

    def getResult(self):
        """
            .   @brief 返回亚像素
            .   @param out: 类型矩阵，大小为 N*12 的 float 矩阵，N为模型点数，【0~1：匹配中心坐标x,y】【2~3：匹配角度和分数】【4~11：矩形四个顶点坐标】
            """
        pass



