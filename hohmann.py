#import ksptools
#import ksptools.orbit as orbit
#import ksptools.util as util

def hohmann_simple(r1, r2):
    from ksptools.util import arcvec, veccos
    from numpy.linalg import norm
    
    dtheta = arcvec(r1, r2)
    cost = veccos(r1,r2)
    r1_len = norm(r1)
    r2_len = norm(r2)
    
    print("dtheta: {}".format(dtheta))
    print("r1: {}".format(r1_len))
    print("r2: {}".format(r2_len))
    print("veccos: {}".format(cost))
    
    rat = r1_len/r2_len
    
    if r1_len > r2_len:
        e = (rat-1)/(cost+rat)
        a = r1_len*(1-e)/(1-e**2)
    elif r1_len < r2_len:
        e = (rat-1)/(cost-rat)
        a = r1_len*(1+e)/(1-e**2)
    return e, a
    
def hohmann_
    
