import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import scipy.interpolate as interpolate
from scipy.interpolate import interp1d

YdataNiCo = np.array ([-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.4])
XdataNiCo = np.array([0.03986 ,0.0904, 0.18840, 0.36659, 0.58265, 0.78241, 0.90809, 0.9842])

XdataNiCr = np.array ([0.023, 0.0645, 0.1139,  0.3555,  0.5117, 0.6343, 0.8140, 0.9247, 0.9443])
YdataNiCr = np.array ([-0.8, -0.7, -0.65, -0.55, -0.5, -0.4, -0.1, 0.3, 0.4])

def Deltamunicr (X, a, b, c, d, e, f, g, h,i):
	return 1500*np.log(X/(1-X)) + a + b*X + c*X*X + d*X*X*X + e*X*X*X*X + f*X*X*X*X*X + g*X*X*X*X*X*X + h*X*X*X*X*X*X*X + i*X*X*X*X*X*X*X*X
	

def DeltaMuNiCr (Ni):
	popt, pcov = curve_fit(Deltamunicr, XdataNiCr, YdataNiCr)
	fig, ax = plt.subplots()
	N=2000
	xmin, xmax = XdataNiCr.min(), XdataNiCr.max()
	xx = np.linspace(xmin, xmax, N)
	t, c, k = interpolate.splrep(XdataNiCr, Deltamunicr(XdataNiCr, *popt), s=0.0001, k=4)
	spline = interpolate.BSpline(t, c, k, extrapolate=False)
	return spline (Ni)
	
#print (DeltaMuNiCr (0.4))


def Deltamunico (X, a, b, c, d, e, f, g, h):
	return 1500*np.log(X/(1-X)) + a + b*X + c*X*X + d*X*X*X + e*X*X*X*X + f*X*X*X*X*X + g*X*X*X*X*X*X + h*X*X*X*X*X*X*X
	
def DeltaMuNiCo (Ni):
	popt, pcov = curve_fit(Deltamunico, XdataNiCo, YdataNiCo)
	N=2000
	xmin, xmax = XdataNiCo.min(), XdataNiCo.max()
	xx = np.linspace(xmin, xmax, N)
	t, c, k = interpolate.splrep(XdataNiCo, Deltamunico(XdataNiCo, *popt), s=0.00001, k=3)
	spline = interpolate.BSpline(t, c, k, extrapolate=False)
	return spline (Ni)
	
#print (DeltaMuNiC0 (0.4))
