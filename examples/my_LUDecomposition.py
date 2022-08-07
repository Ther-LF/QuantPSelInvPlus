import numpy as np
import pandas as pd
 
np.random.seed(2)
def LU_decomposition(A):
    n=len(A[0])
    L = np.zeros([n,n])
    U = np.zeros([n, n])
    for i in range(n):
        L[i][i]=1
        if i==0:
            U[0][0] = A[0][0]
            for j in range(1,n):
                U[0][j]=A[0][j]
                L[j][0]=A[j][0]/U[0][0]
        else:
                for j in range(i, n):#U
                    temp=0
                    for k in range(0, i):
                        temp = temp+L[i][k] * U[k][j]
                    U[i][j]=A[i][j]-temp
                for j in range(i+1, n):#L
                    temp = 0
                    for k in range(0, i ):
                        temp = temp + L[j][k] * U[k][i]
                    L[j][i] = (A[j][i] - temp)/U[i][i]
        print(i)
    return L,U
 
if __name__ == '__main__': 
    #A=np.random.randint(1,10,size=[3,3])  #注意A的顺序主子式大于零
    A=np.loadtxt(open("H180.csv","rb"),delimiter=",")   #举一个例子
    L,U=LU_decomposition(A)
    df = pd.DataFrame(L)
    df.to_csv('L180.csv')
