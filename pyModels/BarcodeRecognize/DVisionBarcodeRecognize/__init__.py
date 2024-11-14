from DVisionBarcodeRecognize import DVisionBarcodeRecognize
from enum import Enum

# 创建模型
class BarcodeRecognize():
    def __new__(self):
        return DVisionBarcodeRecognize.BarcodeRecognize()

    def setFormats(self, formats):
        """
            .   @brief 设置要解码的格式。 如果“无”，则解码所有格式。
            .   @param formats: [解码的格式] [BarcodeFormats] [默认值：BarcodeFormats{}]
            """
        pass
    def setTextMode(self, text_mode):
        """
            .   @brief 指定控制原始字节内容如何转码为结果中的文本的 TextMode
            .   @param text_mode: [类型：TextMode] [默认值：TextMode.HRI]
            """
        pass
    def setBinarizer(self, binarizer):
        """
            .   @brief 二值化器用于在解码条形码之前转换图像
            .   @param binarizer: [类型：Binarizer] [默认值：Binarizer.LocalAverage]
            """
        pass
    def setEanAddOnSymbol(self, ean_add_on_symbol):
        """
            .   @brief  指定在扫描 EAN/UPC 代码时是忽略、读取还是需要 EAN-2/5 附加符号。
            .   @param autoThreshold:  [类型：EanAddOnSymbol] [默认值：EanAddOnSymbol.Ignore]
            """
        pass
    def setRotate(self, try_rotate):
        """
            .   @brief True:解码器在任何方向搜索条形码; False:不会搜索 90° / 270° 旋转的条形码。
            """
        pass
    def setDownscale(self, try_downscale):
        """
            .   @brief True:解码器会扫描输入的缩小版本；False:只会在提供的分辨率中搜索
            """
        pass
    def setPure(self, is_pure):
        """
            .   @brief 如果输入仅包含完美对齐的条形码（生成的图像），则设置为 True,在这种情况下加快检测速度。
            """
        pass
    def recognize(self, _image):
        """
            .   @brief 检测单个码
            .   @param _image: 图像矩阵
            """
        pass
    def getResult(self):
        """
            .   @brief 返回结果，单个码结果
            .   @param out: [类型：Result]
            """
        pass

    def recognizes(self, _image):
        """
            .   @brief 检测多个个码
            .   @param _image: 图像矩阵
            """
        pass

    def getResults(self):
        """
            .   @brief 返回结果，多个码结果
            .   @param out: [类型：Results]
            """
        pass

# Enum BarcodeFormat
m_BarcodeFormat_NONE =DVisionBarcodeRecognize.BarcodeFormat.NONE
m_BarcodeFormat_Aztec =DVisionBarcodeRecognize.BarcodeFormat.Aztec
m_BarcodeFormat_Codabar =DVisionBarcodeRecognize.BarcodeFormat.Codabar
m_BarcodeFormat_Code39 =DVisionBarcodeRecognize.BarcodeFormat.Code39
m_BarcodeFormat_Code93 =DVisionBarcodeRecognize.BarcodeFormat.Code93
m_BarcodeFormat_Code128 =DVisionBarcodeRecognize.BarcodeFormat.Code128
m_BarcodeFormat_DataMatrix =DVisionBarcodeRecognize.BarcodeFormat.DataMatrix
m_BarcodeFormat_EAN8 =DVisionBarcodeRecognize.BarcodeFormat.EAN8
m_BarcodeFormat_EAN13 =DVisionBarcodeRecognize.BarcodeFormat.EAN13
m_BarcodeFormat_ITF =DVisionBarcodeRecognize.BarcodeFormat.ITF
m_BarcodeFormat_MaxiCode =DVisionBarcodeRecognize.BarcodeFormat.MaxiCode
m_BarcodeFormat_PDF417 =DVisionBarcodeRecognize.BarcodeFormat.PDF417
m_BarcodeFormat_QRCode =DVisionBarcodeRecognize.BarcodeFormat.QRCode
m_BarcodeFormat_MicroQRCode =DVisionBarcodeRecognize.BarcodeFormat.MicroQRCode
m_BarcodeFormat_DataBar =DVisionBarcodeRecognize.BarcodeFormat.DataBar
m_BarcodeFormat_DataBarExpanded =DVisionBarcodeRecognize.BarcodeFormat.DataBarExpanded
m_BarcodeFormat_UPCA =DVisionBarcodeRecognize.BarcodeFormat.UPCA
m_BarcodeFormat_UPCE =DVisionBarcodeRecognize.BarcodeFormat.UPCE
m_BarcodeFormat_LinearCodes =DVisionBarcodeRecognize.BarcodeFormat.LinearCodes
m_BarcodeFormat_MatrixCodes =DVisionBarcodeRecognize.BarcodeFormat.MatrixCodes

# Enum Binarizer
m_Binarizer_BoolCast =DVisionBarcodeRecognize.Binarizer.BoolCast
m_Binarizer_LocalAverage =DVisionBarcodeRecognize.Binarizer.LocalAverage
m_Binarizer_FixedThreshold =DVisionBarcodeRecognize.Binarizer.FixedThreshold
m_Binarizer_GlobalHistogram =DVisionBarcodeRecognize.Binarizer.GlobalHistogram

# Enum EanAddOnSymbol
m_EanAddOnSymbol_Ignore =DVisionBarcodeRecognize.EanAddOnSymbol.Ignore
m_EanAddOnSymbol_Require =DVisionBarcodeRecognize.EanAddOnSymbol.Require
m_EanAddOnSymbol_Read =DVisionBarcodeRecognize.EanAddOnSymbol.Read

# Enum ContentType
m_ContentType_Text =DVisionBarcodeRecognize.ContentType.Text
m_ContentType_Binary =DVisionBarcodeRecognize.ContentType.Binary
m_ContentType_Mixed =DVisionBarcodeRecognize.ContentType.Mixed
m_ContentType_GS1 =DVisionBarcodeRecognize.ContentType.GS1
m_ContentType_ISO15434 =DVisionBarcodeRecognize.ContentType.ISO15434
m_ContentType_UnknownECI =DVisionBarcodeRecognize.ContentType.UnknownECI

# Enum TextMode
m_TextMode_Plain =DVisionBarcodeRecognize.TextMode.Plain
m_TextMode_ECI =DVisionBarcodeRecognize.TextMode.ECI
m_TextMode_HRI =DVisionBarcodeRecognize.TextMode.HRI
m_TextMode_Hex =DVisionBarcodeRecognize.TextMode.Hex
m_TextMode_Escaped =DVisionBarcodeRecognize.TextMode.Escaped



