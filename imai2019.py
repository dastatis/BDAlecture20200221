
import pandas as pd
import numpy as np
import math

X = pd.read_csv('C:/Users/dhojo/Desktop/X.csv')

n = len(X)
y = pd.read_csv('C:/Users/dhojo/Desktop/y.csv')

r_0 = 0.06
s_0 = 6
Q_0 = np.diag(np.array([r_0,s_0]))
a_0 = 6
b_0 = 600**2

## the exact value of marginal likelihood

M = np.dot(X.T,X) + Q_0
R = np.identity(n) - np.dot(np.dot(X,np.linalg.inv(M)),X.T)
ML = math.pi**(-n/2) * b_0**(a_0/2) * (math.gamma((n+a_0)/2) / math.gamma(a_0/2)) * \
     (np.linalg.det(Q_0)**(1/2) / np.linalg.det(M)**(1/2)) * (np.dot(np.dot(y.T,R),y)+b_0)**(-(n+a_0)/2)
math.log(ML)

logML = -(n/2)*math.log(math.pi) + (a_0/2)*math.log(b_0) + math.log( math.gamma((n+a_0)/2)/math.gamma(a_0/2) ) + \
  (1/2)*math.log( np.linalg.det(Q_0)/np.linalg.det(M) ) - ((n+a_0)/2)*math.log(np.dot(np.dot(y.T,R),y)+b_0)
logML


