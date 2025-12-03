import numpy as np


def CreateConfig(Nb, Distance, type):
    # Distance = diamètre + d (input)
    coord = np.zeros((Nb, 2))
    if type=="Nb2_X" : 
        coord[0, :] = [0, 0]
        coord[1, :] = [Distance, 0]
    elif type=="Nb2_Y" : 
        coord[0, :] = [0, 0]
        coord[1, :] = [0, Distance]
    elif type=="Nb3_X" : 
        coord[0, :] = [0, 0]
        coord[1, :] = [Distance, 0]
        coord[2, :] = [-(Distance), 0]
    elif type=="Nb3_Y" : 
        coord[0, :] = [0, 0]
        coord[1, :] = [0, Distance]
        coord[2, :] = [0, -(Distance)]
    elif type=="Nb3_T" : 
        coord[0, :] = [(Distance)*np.sqrt(3)/2, 0]
        coord[1, :] = [0, (Distance)/2]
        coord[2, :] = [0, -(Distance)/2]
    elif type=="Nb3_Tinv" : 
        coord[0, :] = [0, 0]
        coord[1, :] = [(Distance)*np.sqrt(3)/2, (Distance)/2]
        coord[2, :] = [(Distance)*np.sqrt(3)/2, -(Distance)/2]
    elif type=="Nb4_C" : 
        coord[0, :] = [-(Distance/2), -(Distance/2)]
        coord[1, :] = [-(Distance/2), (Distance/2)]
        coord[2, :] = [(Distance/2), (Distance/2)]
        coord[3, :] = [(Distance/2), -(Distance/2)]
    elif type=="Nb9_C" : 
        coord[0, :] = [0,0]
        coord[1, :] = [Distance, 0]
        coord[2, :] = [Distance, -Distance]
        coord[3, :] = [0, -Distance]
        coord[4, :] = [-Distance, -Distance]
        coord[5, :] = [-Distance, 0]
        coord[6, :] = [-Distance, Distance]
        coord[7, :] = [0, Distance]
        coord[8, :] = [Distance, Distance]
    return coord[:,:]
