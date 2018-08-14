"""
errors in:
Gaussian fit
Lorentzian fit
Voigt fit

for collapsed x and y coords
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit,leastsq
from astropy.modeling import models,fitting
import lmfit

# def gauss(x,a,mu,sigma):
#     """ Return Gaussian line shape at x with HWHM alpha """
#     
#     return (np.e**(-(x-mu)**2/(2*sigma*sigma)))/(np.sqrt(2*np.pi*sigma*sigma))

def FWHM(s,optim,model):
    
    y = model(np.where(s==np.max(s))[0][0],*optim)/2
    
    first,sec = None,None
    
    width = []
    
    for i in range(len(s)):
        
        if(i):
            if (sec is not None):
                if(first<y<sec)or(first>y>sec):
                    width.append(i)
            
            sec = first.copy()
        
        first = model(i,*optim)
    
    return width
    

def gauss(x,a,mu,sigma):
    return (a/np.sqrt(2*np.pi*sigma*sigma))*np.exp(-(x-mu)**2/(2*sigma**2))

def lorentz(x,I,x0,gamma):
    return I * gamma**2 / ((x-x0)**2 + gamma**2)

def voigt(x,amplitude_L,x_0,fwhm_L,fwhm_G):
    
    A = np.array([-1.2150, -1.3509, -1.2150, -1.3509])
    B = np.array([1.2359, 0.3786, -1.2359, -0.3786])
    C = np.array([-0.3085, 0.5906, -0.3085, 0.5906])
    D = np.array([0.0210, -1.1858, -0.0210, 1.1858])
    
    sqrt_ln2 = np.sqrt(np.log(2))
    X = (x - x_0) * 2 * sqrt_ln2 / fwhm_G
    X = np.atleast_1d(X)[..., np.newaxis]
    Y = fwhm_L * sqrt_ln2 / fwhm_G
    Y = np.atleast_1d(Y)[..., np.newaxis]

    V = np.sum((C * (Y - A) + D * (X - B))/(((Y - A) ** 2 + (X - B) ** 2)), axis=-1)

    return (fwhm_L * amplitude_L * np.sqrt(np.pi) * sqrt_ln2 / fwhm_G) * V

# arr = np.random.normal(100,3,(100,100))
fiz = fits.open('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/HK180725_041_img.fits')

arr = fiz[0].data
fiz.close()

##Gaussian

r,c = arr.shape

x = np.array(np.hsplit(arr,c))
y = np.array(np.vsplit(arr,r))

xra = range(0,c)
yra = range(0,r)

sum_x = np.array([np.sum(i) for i in x])
sum_y = np.array([np.sum(j) for j in y])

mux  = np.mean(sum_x)
sigx = sum_x.std()

optim_gx, covar_gx = curve_fit(gauss,xra,sum_x,p0=[2,mux,sigx])

gx_init = models.Gaussian1D(amplitude=100, mean=mux, stddev=sigx)
fit_gx = fitting.LevMarLSQFitter()
gx = fit_gx(gx_init, xra, sum_x)

print('='*50)
print('Gaussian')
print('='*50)

plt.figure()
plt.plot(xra,sum_x,'r*',label='data')
plt.plot(xra,gauss(xra,*optim_gx),'b-',label='scipy.curve_fit')
plt.plot(xra,gx(xra),'g--',label='astropy.fitting')
plt.title('Compressed X Gaussian Fit')
plt.xlabel('Position')
plt.ylabel('Events')
plt.legend()

print('x-direction')
print('-'*20)
print('Amplitude:',optim_gx[0])
print('Mean:     ',optim_gx[1])
print('Sigma:    ',optim_gx[2])
print('Error:    ',np.sqrt(np.diag(covar_gx)))
print('Peak:     ',gauss(np.where(sum_x==np.max(sum_x))[0][0],*optim_gx))
print('FWHM:     ',FWHM(sum_x,optim_gx,gauss))

muy  = np.mean(sum_y)
sigy = sum_y.std()

optim_gy, covar_gy = curve_fit(gauss,yra,sum_y,p0=[4,muy,sigy])

gy_init = models.Gaussian1D(amplitude=120, mean=muy, stddev=sigy)
fit_gy = fitting.LevMarLSQFitter()
gy = fit_gy(gy_init, yra, sum_y)

plt.figure()
plt.plot(yra,sum_y,'r*',label='data')
plt.plot(yra,gauss(yra,*optim_gy),'b-',label='scipy.curve_fit')
plt.plot(yra,gy(yra),'g--',label='astropy.fitting')
plt.xlabel('Position')
plt.ylabel('Events')
plt.title('Compressed Y Gaussian Fit')
plt.legend()

print('\ny-direction')
print('-'*20)
print('Amplitude:',optim_gy[0])
print('Mean:     ',optim_gy[1])
print('Sigma:    ',optim_gy[2])
print('Error:    ',np.sqrt(np.diag(covar_gy)))
print('Peak:     ',gauss(np.where(sum_y==np.max(sum_y))[0][0],*optim_gy))
print('FWHM:     ',FWHM(sum_y,optim_gy,gauss))

##Lorentzian

'''max values'''
maxx = np.max(sum_x)
maxy = np.max(sum_y)
'''max locations'''
max_locx = np.where(sum_x == maxx)[0]
max_locy = np.where(sum_y == maxy)[0]
'''half maxes'''
hmaxx = 0.5 * maxx
hmaxy = 0.5 * maxy
'''some BrenDan math'''
v_x = optim_gx[2]*np.sqrt(-2*np.log(hmaxx/optim_gx[0]))+optim_gx[1]
v_y = optim_gy[2]*np.sqrt(-2*np.log(hmaxy/optim_gy[0]))+optim_gy[1]
'''half-width half maxes'''
hwhmx = abs(max_locx - v_x)
hwhmy = abs(max_locy - v_y)
'''fit the lorentzian with SciPy'''
optim_lx, covar_lx = curve_fit(lorentz,xra,sum_x,p0=[1,max_locx,hwhmx])
optim_ly, covar_ly = curve_fit(lorentz,yra,sum_y,p0=[1,max_locy,hwhmy])

'''FWHM and HWHM my way'''
gwx = FWHM(sum_x,optim_gx,gauss)
dif_gx = gwx[1]-gwx[0]
gwy = FWHM(sum_y,optim_gy,gauss)
dif_gy = gwy[1]-gwy[0]

lx_init = models.Lorentz1D(amplitude=120, x_0=max_locx, fwhm=dif_gx)
fit_lx = fitting.LevMarLSQFitter()
lx = fit_lx(lx_init, xra, sum_x)

ly_init = models.Lorentz1D(amplitude=120, x_0=max_locy, fwhm=dif_gy)
fit_ly = fitting.LevMarLSQFitter()
ly = fit_ly(ly_init, yra, sum_y)

print('='*50)
print('Lorentzian')
print('='*50)

plt.figure()
plt.plot(xra,sum_x,'r*',label='data')
plt.plot(xra,lorentz(xra,*optim_lx),'b-',label='scipy.curve_fit')
plt.plot(xra,lx(xra),'g--',label='astropy.fitting')
plt.title('Compressed X Lorentzian Fit')
plt.xlabel('Position')
plt.ylabel('Events')
plt.legend()

print('x-direction')
print('-'*20)
print('Amplitude:',optim_lx[0])
print('Mean:     ',optim_lx[1])
print('Sigma:    ',optim_lx[2])
print('Error:    ',np.sqrt(np.diag(covar_lx)))
print('Peak:     ',lorentz(np.where(sum_x==np.max(sum_x))[0][0],*optim_lx))
print('FWHM:     ',FWHM(sum_x,optim_lx,lorentz))

plt.figure()
plt.plot(yra,sum_y,'r*',label='data')
plt.plot(yra,lorentz(yra,*optim_ly),'b-',label='scipy.curve_fit')
plt.plot(yra,ly(yra),'g--',label='astropy.fitting')
plt.title('Compressed Y Lorentzian Fit')
plt.xlabel('Position')
plt.ylabel('Events')
plt.legend()

print('\ny-direction')
print('-'*20)
print('Amplitude:',optim_ly[0])
print('Mean:     ',optim_ly[1])
print('Sigma:    ',optim_ly[2])
print('Error:    ',np.sqrt(np.diag(covar_ly)))
print('Peak:     ',lorentz(np.where(sum_y==np.max(sum_y))[0][0],*optim_ly))
print('FWHM:     ',FWHM(sum_y,optim_ly,lorentz))

## Voigt

lwx = FWHM(sum_x,optim_lx,lorentz)
dif_lx = lwx[1]-lwx[0]
lwy = FWHM(sum_y,optim_ly,lorentz)
dif_ly = lwy[1]-lwy[0]

vx_init = models.Voigt1D(x_0 = max_locx, amplitude_L = optim_lx[0], \
                         fwhm_L = dif_lx, fwhm_G = dif_gx)
fit_vx = fitting.LevMarLSQFitter()
vx = fit_vx(vx_init, xra, sum_x)

vy_init = models.Voigt1D(x_0 = max_locy, amplitude_L = optim_ly[0], \
                         fwhm_L = dif_ly, fwhm_G = dif_gy)
fit_vy = fitting.LevMarLSQFitter()
vy = fit_vy(vy_init, yra, sum_y)


optim_vx, covar_vx = curve_fit(voigt,xra,sum_x,p0=[1,max_locx,dif_lx,dif_gx])
optim_vy, covar_vy = curve_fit(voigt,yra,sum_y,p0=[1,max_locy,dif_ly,dif_gy])

print('='*50)
print('Voigt')
print('='*50)

print('x-direction')
print('-'*20)
print('Amplitude:',optim_vx[0])
print('Mean:     ',optim_vx[1])
print('Sigma:    ',optim_vx[2])
print('Error:    ',np.sqrt(np.diag(covar_vx)))
print('Peak:     ',voigt(np.where(sum_x==np.max(sum_x))[0][0],*optim_vx))
print('FWHM:     ',FWHM(sum_x,optim_vx,voigt))

print('\ny-direction')
print('-'*20)
print('Amplitude:',optim_vy[0])
print('Mean:     ',optim_vy[1])
print('Sigma:    ',optim_vy[2])
print('Error:    ',np.sqrt(np.diag(covar_vy)))
print('Peak:     ',voigt(np.where(sum_y==np.max(sum_y))[0][0],*optim_vy))
print('FWHM:     ',FWHM(sum_y,optim_vy,voigt))

plt.figure()
plt.plot(xra,sum_x,'r*',label='data')
plt.plot(xra,voigt(xra,*optim_vx),'b-',label='scipy.curve_fit')
plt.plot(xra,vx(xra),'g--',label='astropy.fitting')
plt.title('Compressed X Voigt Fit')
plt.xlabel('Position')
plt.ylabel('Events')
plt.legend()


plt.figure()
plt.plot(yra,sum_y,'r*',label='data')
plt.plot(yra,voigt(yra,*optim_vy),'b-',label='scipy.curve_fit')
plt.plot(yra,vy(yra),'g--',label='astropy.fitting')
plt.title('Compressed Y Voigt Fit')
plt.xlabel('Position')
plt.ylabel('Events')
plt.legend()

plt.show()

