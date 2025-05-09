
from scipy.special import jv, yv
import numpy as np
PI=4.*np.arctan(1.)
N=4
print("N=", N)
X=1
H_1=jv(N, X)+1j*yv(N, X)
H_2=jv(N, X)-1j*yv(N, X)

# print("H_1", H_1)
print("H_2", H_2)

# Hankel_2 = -1j * jv(N, X) * (np.exp(PI*N*1j) -(-1)**N)/ np.sin(PI * N)
Hankel_2 = -1j * (jv(N, X) * np.exp(PI*N*1j)-jv(-N, X))/ np.sin(PI * N)
print("Hankel", Hankel_2)

# print("test", yv(N,X), (-1)**(N+1/2)*jv(-N,X))