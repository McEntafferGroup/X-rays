import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import norm

fitsfile = fits.open('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/HK180725_041_img.fits')

data = np.array(fitsfile[0].data)

r, c = data.shape

x_col = np.array(np.hsplit(data, c))
y_col = np.array(np.vsplit(data, r))

x_sum = np.array([np.sum(i) for i in x_col])
y_sum = np.array([np.sum(j) for j in y_col])

xra = range(0,c)
yra = range(0,r)



mux = np.mean(x_sum)
sigmax = np.std(x_sum)

muy = np.mean(y_sum)
sigmay = np.std(y_sum)

def norm(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2 / (2*sigma**2))
    
def lorentz(x,I,x0,gamma):
    return I * gamma**2 / ((x-x0)**2 + gamma**2)

maxx = np.max(x_sum)
maxy = np.max(y_sum)

max_locx = np.where(x_sum == maxx)[0]
max_locy = np.where(y_sum == maxy)[0]

hmaxx = 0.5 * maxx
hmaxy = 0.5 * maxy

h_locx = []
h_locy = []


optimx, covarx = curve_fit(norm, xra, x_sum, p0=[2,mux,sigmax])
optimy, covary = curve_fit(norm, yra, y_sum, p0=[3,muy,sigmay])
 
v_x = optimx[2] * np.sqrt(-2*np.log(hmaxx/optimx[0]))+optimx[1]
v_y = optimy[2] * np.sqrt(-2*np.log(hmaxy/optimy[0]))+optimy[1]

hwhmx = abs(max_locx - v_x)
hwhmy = abs(max_locy - v_y)

print(hwhmx, hwhmy)

optimxl, covarxl = curve_fit(lorentz, xra, x_sum, p0=[1,max_locx,hwhmx])
optimyl, covaryl = curve_fit(lorentz, yra, y_sum, p0=[1,max_locy,hwhmy])

plt.plot(xra, x_sum, 'r*')
plt.plot(xra, norm(xra, *optimx), 'b-')
plt.title('Gaussian Fit to x')
plt.figure()

plt.plot(yra, y_sum, 'r*')
plt.plot(yra, norm(yra, *optimy), 'b-')
plt.title('Gaussian Fit to y')
plt.figure()

plt.plot(xra, x_sum, 'r*')
plt.plot(xra, lorentz(xra, *optimxl), 'g-')
plt.title('Lorentzian Fit to x')
plt.figure()

plt.plot(yra, y_sum, 'r*')
plt.plot(yra, lorentz(yra, *optimyl), 'g-')
plt.title('lorentzian Fit to y')
plt.show()
 
print(optimx)
print(np.diag(covarx))
print(optimy)
print(np.diag(covary))

