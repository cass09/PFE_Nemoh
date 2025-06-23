import numpy as np

Nfinal=2100
cyl=False


if cyl :
    Ninit=361
    Nr=5
    Ntheta=19
    Nz=14


    facteur=np.sqrt(Nfinal/Ninit)
    Nr=int(Nr*facteur)
    Ntheta=int(Ntheta*facteur)
    Nz=int(Nz*facteur)

    print("Nr, Ntheta, Nz = ", Nr, Ntheta, Nz)
    print("Npanels = ", Ntheta*(Nr+Nz))

else :
    N=13
    Ninit=403


    facteur=np.sqrt(Nfinal/Ninit)
    N=int(N*facteur)

    print("N = ", N)
    print("Npanels = ", 5*N**2)
    print("Npanels 1/2= ", N**2+3*N*int(N/2))
