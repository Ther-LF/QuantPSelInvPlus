from turtle import shape
import pandas as pd 
import numpy as np
from functools import cmp_to_key

def cmp(a, b):
    if abs(a.avg) < abs(b.avg):
        return -1
    else:
        return 1

class node: 
    def __init__(self, row, col, avg):
        self.row = row
        self.col = col
        self.avg = avg



data = pd.read_csv(r'/home/luofan/gitRepo/QuantPSelInvPlus/examples/L180.csv',sep=',')
# print(data.shape)
arr = np.array(data)
superidx = [0, 60, 100, 200, 260, 300, 400, 500, 600, 660, 700, 800, 900, 1000, 1060, 1100, 1200, 1300, 1400, 1600, 1856, 2112, 2368, 2624, 2880, 3136, 3392, 3600]
nodes = []
for r in range(27):
    for c in range(r + 1):
        sum = 0.0
        for i in range(superidx[r], superidx[r + 1]):
            for j in range(superidx[c], superidx[c + 1]):
                sum = sum + arr[i][j]
        avg = sum / (superidx[r + 1] - superidx[r]) / (superidx[c + 1] - superidx[c])
        nodes.append(node(r, c, avg))

nodes.sort(key=cmp_to_key(cmp))
for n in nodes:
    print(n.row, n.col, n.avg)
        


