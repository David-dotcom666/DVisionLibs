import cv2
import DVisionBarcodeRecognize
import numpy as np
# 获取图像

img = cv2.imread("image/bar2.bmp")
img2 = cv2.imread("image/bar4.bmp")
t1 =DVisionBarcodeRecognize.BarcodeRecognize()
t1.setFormats(DVisionBarcodeRecognize.m_BarcodeFormat_NONE)
t1.setBinarizer(DVisionBarcodeRecognize.m_Binarizer_LocalAverage)

# 单个目标
t1.recognize(img)
result2 = t1.getResult()
print("Found barcode:\n Text:    '{}'\n Format:   {}\n Content:  {}\n Position: {}"
		.format(result2.text, result2.format, result2.content_type, result2.position))

print("***********************")

# 多个目标
t1.recognizes(img)
results = t1.getResults()
# 定义线的颜色 (B, G, R)，这里为红色
color = (0, 255, 0)
# 定义线的宽度
thickness = 2
for result in results:
	print("Found barcode:\n Text:    '{}'\n Format:   {}\n Content:  {}\n Position: {}"
		.format(result.text, result.format, result.content_type, result.position))
	position = result.position
	cv2.line(img, (position.top_left.x, position.top_left.y), (position.top_right.x, position.top_right.y), color,
			 thickness)
	cv2.line(img, (position.top_right.x, position.top_right.y), (position.bottom_right.x, position.bottom_right.y), color,
			 thickness)
	cv2.line(img, (position.bottom_right.x, position.bottom_right.y), (position.bottom_left.x, position.bottom_left.y), color,
			 thickness)
	cv2.line(img, (position.top_left.x, position.top_left.y), (position.bottom_left.x, position.bottom_left.y), color,
			 thickness)
	centerx = round((position.top_right.x + position.top_left.x + position.bottom_right.x + position.bottom_left.x) / 2)
	centery = round((position.top_right.y + position.top_left.y + position.bottom_right.y + position.bottom_left.y) / 2)
	print(centerx, centery, result.text)
	cv2.putText(img, result.text, (position.top_left.x, position.top_left.y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255), 1, cv2.LINE_AA)

cv2.imshow("img", img)
cv2.waitKey()