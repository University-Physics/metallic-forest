import cv2
import numpy as np
re_img1 = cv2.imread('Raunkiaer.jpg')
b, g, r = cv2.split(re_img1)
print(b)
