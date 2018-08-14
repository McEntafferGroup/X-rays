import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy
from scipy.optimize import curve_fit,leastsq
from astropy.modeling import models,fitting
import lmfit
from lmfit import minimize, Parameters, Model

def FWHM(s,optim,model):
    
    y = model.eval(params=optim,x=np.where(s==np.max(s))[0][0])/2
    
    first,sec = None,None
    
    width = []
    
    for i in range(len(s)):
        
        if(i):
            if (sec is not None):
                if(first<y<sec)or(first>y>sec):
                    width.append(i)
            
            sec = first.copy()
        
        first = model.eval(params=optim,x=i)
    
    return width

def gauss(x,sigma,mu,a):
    return (a/np.sqrt(2*np.pi*sigma*sigma))*np.exp(-(x-mu)**2/(2*sigma**2))
    
    
def plot_gaussian(arr, *bounds):
    
    try:
        
        arr = arr[:,bounds[0]:bounds[1]]
        arr = arr[bounds[0]:bounds[1]]
        
        xra = np.arange(bounds[0],bounds[1])
        yra = np.arange(bounds[0],bounds[1])
        
        r,c = arr.shape
    
    except:
        
        r,c = arr.shape
        
        xra = np.arange(0,c)
        yra = np.arange(0,r)
    
    x = np.array(np.hsplit(arr,c))
    y = np.array(np.vsplit(arr,r))
    
    xra = np.arange(0,c)
    yra = np.arange(0,r)
    
    sum_x = np.array([np.sum(i) for i in x])
    sum_y = np.array([np.sum(j) for j in y])
    
    max_locx = np.where(sum_x==sum_x.max())[0]
    max_locy = np.where(sum_y==sum_y.max())[0]
    
    sigx = sum_x.std()
    sigy = sum_y.std()
    
    wx = 1./np.sqrt(sum_x)
    wx[wx==np.inf] = 0
    wy = 1./np.sqrt(sum_y)
    wy[wy==np.inf] = 0
    
    ## Gaussian
    
    gmod = lmfit.models.GaussianModel()
    
    pargx = gmod.guess(sum_x, x=xra)
    gx = gmod.fit(sum_x, pargx, x=xra, weights=wx)
    
    pargy = gmod.guess(sum_y, x=yra)
    gy = gmod.fit(sum_y, pargy, x=yra, weights=wy)
    
    print('='*50)
    print('Gaussian')
    print('='*50)
    
    print('\nx-direction')
    print('-'*20)
    print(gx.fit_report())
    
    print('\ny-direction')
    print('-'*20)
    print(gy.fit_report())
    
    plt.figure()
    plt.plot(xra,sum_x,'r*',label='data points')
    plt.plot(xra,gx.init_fit,'g--',label='initial values')
    plt.plot(xra,gx.best_fit,'b-',label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed X Gaussian Fit')
    
    plt.figure()
    plt.plot(yra,sum_y,'r*',       label='data points')
    plt.plot(yra,gy.init_fit,'g--',label='initial values')
    plt.plot(yra,gy.best_fit,'b-', label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed Y Gaussian Fit')
    
    plt.show()

def plot_lorentzian(arr, *bounds):

    try:
        
        arr = arr[:,bounds[0]:bounds[1]]
        arr = arr[bounds[0]:bounds[1]]
        
        xra = np.arange(bounds[0],bounds[1])
        yra = np.arange(bounds[0],bounds[1])
        
        r,c = arr.shape
    
    except:
        
        r,c = arr.shape
        
        xra = np.arange(0,c)
        yra = np.arange(0,r)
    
    x = np.array(np.hsplit(arr,c))
    y = np.array(np.vsplit(arr,r))
    
    xra = np.arange(0,c)
    yra = np.arange(0,r)
    
    sum_x = np.array([np.sum(i) for i in x])
    sum_y = np.array([np.sum(j) for j in y])
    
    max_locx = np.where(sum_x==sum_x.max())[0]
    max_locy = np.where(sum_y==sum_y.max())[0]
    
    sigx = sum_x.std()
    sigy = sum_y.std()
    
    wx = 1./np.sqrt(sum_x)
    wx[wx==np.inf] = 0
    wy = 1./np.sqrt(sum_y)
    wy[wy==np.inf] = 0
    
    ## Lorentzian
    
    lmod = lmfit.models.LorentzianModel()
    
    parlx = lmod.guess(sum_x, x=xra)
    lx = lmod.fit(sum_x, parlx, x=xra, weights=wx)
    
    parly = lmod.guess(sum_y, x=yra)
    ly = lmod.fit(sum_y, parly, x=yra, weights=wy)
    
    print('='*50)
    print('Lorentzian')
    print('='*50)
    
    print('\nx-direction')
    print('-'*20)
    print(lx.fit_report())
    
    print('\ny-direction')
    print('-'*20)
    print(ly.fit_report())
    
    plt.figure()
    plt.plot(xra,sum_x,'r*',label='data points')
    plt.plot(xra,lx.init_fit,'g--',label='initial values')
    plt.plot(xra,lx.best_fit,'b-',label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed X Lorentzian Fit')
    
    plt.figure()
    plt.plot(yra,sum_y,'r*',       label='data points')
    plt.plot(yra,ly.init_fit,'g--',label='initial values')
    plt.plot(yra,ly.best_fit,'b-', label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed Y Lorentzian Fit')
    
    plt.show()

def plot_voigt(arr, *bounds):
    
    try:
        
        arr = arr[:,bounds[0]:bounds[1]]
        arr = arr[bounds[0]:bounds[1]]
        
        xra = np.arange(bounds[0],bounds[1])
        yra = np.arange(bounds[0],bounds[1])
        
        r,c = arr.shape
    
    except:
        
        r,c = arr.shape
        
        xra = np.arange(0,c)
        yra = np.arange(0,r)
    
    
    x = np.array(np.hsplit(arr,c))
    y = np.array(np.vsplit(arr,r))
    
    sum_x = np.array([np.sum(i) for i in x])
    sum_y = np.array([np.sum(j) for j in y])
    
    max_locx = np.where(sum_x==sum_x.max())[0]
    max_locy = np.where(sum_y==sum_y.max())[0]
    
    sigx = sum_x.std()
    sigy = sum_y.std()
    
    wx = 1./np.sqrt(sum_x)
    wx[wx==np.inf] = 0
    wy = 1./np.sqrt(sum_y)
    wy[wy==np.inf] = 0
    
    ## Voigt
    
    vmod = lmfit.models.VoigtModel()
    
    parvx = vmod.guess(sum_x, x=xra)
    vx = vmod.fit(sum_x, parvx, x=xra, weights=wx)
    
    parvy = vmod.guess(sum_y, x=yra)
    vy = vmod.fit(sum_y, parvy, x=yra, weights=wy)
    
    print('='*50)
    print('Voigt')
    print('='*50)
    
    print('\nx-direction')
    print('-'*20)
    print(vx.fit_report())
    
    print('\ny-direction')
    print('-'*20)
    print(vy.fit_report())
    
    plt.figure()
    plt.plot(xra,sum_x,'r*',label='data points')
    plt.plot(xra,vx.init_fit,'g--',label='initial values')
    plt.plot(xra,vx.best_fit,'b-',label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed X Voigt Fit')
    
    plt.figure()
    plt.plot(yra,sum_y,'r*',       label='data points')
    plt.plot(yra,vy.init_fit,'g--',label='initial values')
    plt.plot(yra,vy.best_fit,'b-', label='best fit values')
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Events')
    plt.title('Compressed Y Voigt Fit')

    plt.show()





fiz = fits.open('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/HK180725_041_img.fits')

arr = fiz[0].data
fiz.close()

plot_gaussian(arr)

plot_lorentzian(arr)

plot_voigt(arr,175,300)
































