import numpy as np
import matplotlib.pyplot as plt
from pytest import mark

# plt.rcParams['figure.figsize'] = (60.0, 30.0)

# plt.rcParams['savefig.dpi'] = 60000 #图片像素
# plt.rcParams['figure.dpi'] = 300 #分辨率
# plt.figure(figsize=(60,30), dpi=3000)
data = np.loadtxt(open("H180.csv","rb"),delimiter=",") 
# x = [0,1,2,3]
# y = [0,1,2,3]
# alpha = [-10,0,-10,-5]
# sz = 20
# sc = plt.scatter(x,y,s = sz,c=alpha,marker='o',cmap='Blues_r')
# plt.colorbar(sc)
# plt.show()
# X = [[0, -1, 2, 3],[1, -2, 3, 4]]
ax = plt.matshow(data, cmap=plt.cm.Greens)
# plt.colorbar(ax.colorbar, fraction=0.025)
plt.title("matrix A")
plt.show()