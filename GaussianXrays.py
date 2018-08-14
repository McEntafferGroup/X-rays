"""
run in Python 3.6 32-bit

step 1: Make image with Gaussian noise
step 2: add x-rays
step 3: find the x-rays
step 4: make x-ray splits
step 5: find the splits
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy import ndimage
from scipy.signal import *
from astropy.io import fits

def make_gaussian(size, peak, sigma, show):
    """ 
    size:  (x,y)
    peak:  counts in adu
    sigma: the sigma of the gaussian noise
    show:  boolean of whether to display the image or not    
    
    retruns: 2D np-array of image pixels
    """
    
    noise = np.random.normal(peak,sigma,size)
    
    if (show):
        fig = plt.figure()
        ax1 = plt.subplot(121)
        plt.imshow(noise, origin='lower')
        plt.colorbar()
        
        ax2 = plt.subplot(122)
        plt.hist(noise.flatten(),bins=300)
        
        plt.show()
    
    return noise
    
def add_xrays(image, num, peak, sigma, show):
    """
    image: the np-array of values with gaussian noise
    num:   number of x-rays to be added
    peak:  counts in adu of peak x-ray events
    sigma: the sigma of the x-ray events
    show:  boolean of whether to display the image or not
    
    retruns: 2D np-array of image pixels
    """
    events = np.random.normal(peak,sigma,num)
    x_positions = np.random.random_integers(0,len(image[0])-1,num)
    y_positions = np.random.random_integers(0,len(image)-1,   num)
    
    data = np.column_stack((x_positions, y_positions, events))
    
    print('Single Events')
    print(tabulate(data, headers=['x-pos', 'y-pos','energy'], tablefmt='fancy_grid'))
    
    events_image = image
    
    for i in range(num):
        events_image[y_positions[i]][x_positions[i]] += events[i]
    
    if (show):
        fig = plt.figure()
        plt.imshow(events_image, origin='lower')
        plt.colorbar()
        plt.show()
    
    return events_image


def find_xray(image, thresh):
    """
    image:  the image in which we will find x-ray events
    thresh: threshold of sigmas above which are events
    
    returns: void (will print a nice table)
    """
    
    peak  = image.mean()
    sigma = image.std()
    
    """things past 5*sigma are certainly events"""
    subtracted = (image - peak) - (thresh*sigma)
    
    events = subtracted[subtracted>0]
    pos = list()
    
    for i in events:
        pos.append(np.where(subtracted == i))
    
    
    """reshaping the array to be useful"""
    pos = np.array(pos)
    pos = np.reshape(pos, (len(events),2))
    pos.T[[0, 1]] = pos.T[[1, 0]]
    
    """reshape events into a column"""
    events = np.reshape(( (events + peak) + (thresh*sigma) ), (len(events),1))
    
    """nice pretty data table"""
    data = np.column_stack((pos,events))
    
    """nice pretty output"""
    print('Found Events')
    print(tabulate(data, headers=['x-pos', 'y-pos','energy'], tablefmt='fancy_grid'))


def add_split_events(image, num, type, peak, sigma, show):
    """
    image: the image to add events to
    num:   number of events to be added
    type:  2,4,5 or all splits (2,4,5,3 respectively)
    peak:  counts in adu of peak x-ray events
    sigma: the sigma of the x-ray events
    show:  boolean of whether to display the image or not
    
    retruns: 2D np-array of image pixels
    """
    
    """handle type first"""
    if (type == 2):
        """x or y axis to be split along, per event"""
        xy = np.random.random_integers(0,1,num)
        
        """make the events"""
        events = np.random.normal(peak,sigma,num) 
        splits = np.random.rand(num)
        
        """coefficients of events to split"""
        half1  = splits
        half2  = 1-half1
        
        """place the events"""
        for i in range(num):
            """make it generalized"""
            dx = (0)if(xy[i])else(1)
            dy = (1)if(xy[i])else(0)
            """generalized code for either value"""
            x_pos = np.random.random_integers(0,len(image[0])-1-dx)
            y_pos = np.random.random_integers(0,len(image)   -1-dy)
            
            """fill first pixel"""
            image[y_pos]   [x_pos]    += half1[i]*events[i]
            
            """fill second"""
            image[y_pos+dy][x_pos+dx] += half2[i]*events[i]
            
            print('\n2 Split Event')
            print('x:',x_pos)
            print('y:',y_pos)
            print('E:',events[i])
    
    
    if (type == 4):
        """make the events"""
        events = np.random.normal(peak/4,sigma,num) 
        """coefficients of events to split"""
        bit1 = np.random.normal(1,0.05,num) 
        bit2 = np.random.normal(1,0.05,num) 
        bit3 = np.random.normal(1,0.05,num) 
        bit4 = np.random.normal(1,0.05,num) 
        
        x_pos = np.random.random_integers(0,len(image[0])-2,num)
        y_pos = np.random.random_integers(0,len(image)   -2,num)
        
        for i in range(num):
            """fill first pixel"""
            image[y_pos[i]]  [x_pos[i]]   += bit1[i]*events[i]
            """fill second"""
            image[y_pos[i]+1][x_pos[i]]   += bit2[i]*events[i]
            """fill third"""
            image[y_pos[i]]  [x_pos[i]+1] += bit3[i]*events[i]
            """fill fourth"""
            image[y_pos[i]+1][x_pos[i]+1] += bit4[i]*events[i]
            
            print('\n4 Split Event')
            print('x:',x_pos[i])
            print('y:',y_pos[i])
            print('E:',events[i]*4)
    
    
    if (type == 5):
        """make events"""
        event = np.random.normal(peak/7*3,sigma,num) 
        arms  = np.random.normal(peak/7,  sigma,num) 
        
        x_pos = np.random.random_integers(1,len(image[0])-2,num)
        y_pos = np.random.random_integers(1,len(image)   -2,num)
        
        for i in range(num):
            """fill main pixel"""
            image[y_pos[i]]  [x_pos[i]]   += event[i]
            
            """fill upper"""
            image[y_pos[i]-1][x_pos[i]]   += arms[i]*np.random.normal(1,0.05)
            
            """fill left"""
            image[y_pos[i]]  [x_pos[i]-1] += arms[i]*np.random.normal(1,0.05)
            
            """fill lower"""
            image[y_pos[i]+1][x_pos[i]]   += arms[i]*np.random.normal(1,0.05)
            
            """fill right"""
            image[y_pos[i]]  [x_pos[i]+1] += arms[i]*np.random.normal(1,0.05)
            
            print('\nCross Split Event')
            print('x:',x_pos[i])
            print('y:',y_pos[i])
            print('E:',event[i]*7/3)
    
    if (type == 3):
        """MAKE IT RECURSIVE, WOO!!!"""
        for i in range(num):
            choice = np.random.random_integers(0,2)
            
            if (choice == 1):
                """call this function with a 2"""
                image = add_split_events(image, 1, 2, peak, sigma, 0)
                
            elif (choice == 2):
                """call this function with a 5"""
                image = add_split_events(image, 1, 5, peak, sigma, 0)
            
            else:
                """call this function with a 4"""
                image = add_split_events(image, 1, 4, peak, sigma, 0)
                
        
    if (show):
        fig = plt.figure()
        plt.imshow(image, origin='lower')
        plt.colorbar()
        plt.show()
    
    return image


def find_all_events(image,thresh):
    """
    image:  the image in which we will find x-ray events
    thresh: threshold of sigmas above which are events
    
    returns: image with consolidated events
    """
    image = event_consolidator(image,thresh)
    
    find_xray(image,thresh)
    
    return image

image = make_gaussian((200,200),100,5,0)

image = add_xrays(image, 0, 1000, 10, 0)
    
image = add_split_events(image, 5000, 3, 1000, 10, 1)

# im = fits.open('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/image.fits')
# 
# image = im[0].data
# 
# im.close()

image = find_all_events(image,3)

fig = plt.figure()
plt.imshow(image, origin='lower')
plt.colorbar()
plt.show()