from numpy import *
from scipy import io


def FindCentroidLaserSpotUsingFourier(im):
    """
    function finds centroid of a white laserspot within a dark background
    using fourier algorithm in a two dimensional grayscale image like the one
    in the file 'shot.mat', that you can import and then try the function
    with script test.m. the script loads an image, calculates centroid, shows
    image and position of the centroid with a green haircross.
    """
    
    #image to a double
    vec = double(im)
    #shape of the image
    rows, cols = shape(im)
    
    #1D column-vector of length rows
    i = arange(rows)
    SIN_A = sin(i * 2 * pi / (rows-1))
    COS_A = cos(i * 2 * pi / (rows-1))
    
    #1D vector of length columns
    j = arange(cols)
    j = j[:,newaxis]
    SIN_B = sin(j * 2 * pi / (cols-1))
    COS_B = cos(j * 2 * pi / (cols-1))
    
    a = sum(matmul(COS_A, vec))
    b = sum(matmul(SIN_A, vec))
    c = sum(matmul(vec, COS_B))
    d = sum(matmul(vec, SIN_B))
    
    if (a > 0):
        if (b > 0):
            rphi = 0
        else:
            rphi = 2 * pi
    
    else:
        rphi = pi
    
    
    
    if (c > 0):
        if (d > 0):
            cphi = 0
        else:
            cphi = 2 * pi
    
    else:
        cphi = pi
    
    y = ((b/a) + rphi) * (rows-1) / 2 / pi 
    x = ((d/c) + cphi) * (cols-1) / 2 / pi 
    
    
    return (y,x)
    


if __name__ == '__main__':
    
    abc = io.loadmat('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/Test Images/Test Images/gaussPixelTest')
    
    abc = abc['gFilter']
    
    y,x = FindCentroidLaserSpotUsingFourier(abc)
    
    print('y =',y)
    print('x =',x)















