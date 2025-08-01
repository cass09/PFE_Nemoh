import numpy as np


# k=np.linspace(0, 1, 20)
# k=2*np.pi/(30*3)
a=3
L_wave=3*a
#L_wave=15.3
k=2*np.pi/L_wave
#k=1
g=9.81
#h=50*a/3
h=30
w=np.sqrt(g*k*np.tanh(k*h))
print("h=", h)
print("k=", k, "lambda=", L_wave)
print("w=", w)
print("f=", w/(2*np.pi), "T=", (2*np.pi)/w)
# print("w in ", min(w), max(w))
# print("f in ", min(w)/(2*np.pi), max(w)/(2*np.pi))
