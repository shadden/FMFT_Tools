import matplotlib.pyplot as plt
from ctypes import *
import numpy as np

#LIBPATH = "/Users/samuelhadden/19_FMFT"
LIBPATH = "/projects/b1002/shadden/12_FMFT_tools"

def get_ctype_ptr(dtype,dim,**kwargs):
	return np.ctypeslib.ndpointer(dtype=dtype,ndim=dim,flags='C',**kwargs)
	
p1d=get_ctype_ptr(np.float,1)
p2d = get_ctype_ptr(np.float,2)
pd = POINTER(c_double)
ppd = POINTER(pd)
p1dInt = get_ctype_ptr(np.int,1)

def nearest_pow2(x):
	return int(2**np.floor(np.log2(x)))
	
class libwrapper(object):
	""" A wrapper class for the FMFT library.  Gives access to the function fmft"""
	def __init__(self):
		
		self.lib = CDLL("%s/libfmft.so" % LIBPATH)
		
		self._fmft = self.lib.fmftWrap
		self._fmft.argtypes =[p2d, c_int, c_double, c_double, c_int, p2d, c_long]
		self._fmft.restype = c_int
		def check_errors(ret, func, args):
			if ret<=0:
				raise RuntimeError("FMFT returned error code %d for the given arguments"%ret)
			return ret
		self._fmft.errcheck = check_errors
	
	
		self._array2D = self.lib.array2D
		self._array2D.argtypes = [p2d, c_int]
		self._array2D.restype = c_double

	def example(self,data):
		n = len(data)
		return self._array2D(data,n)
		
	def fmft_full(self, nfreq, minfreq, maxfreq, flag, inpt, ndata):
		"""
		In the output array **output: output[3*flag-2][i], output[3*flag-1][i] 
		and output[3*flag][i] are the i-th frequency, amplitude and phase; nfreq is the 
		number of frequencies to be computed (the units are rad/sep, where sep is the 
		`time' separation between i and i+1. The algorithm is  

		Basic Fourier Transform algorithm           if   flag = 0;   not implemented   
		Modified Fourier Transform                  if   flag = 1;
		Frequency Modified Fourier Transform        if   flag = 2;
		FMFT with additional non-linear correction  if   flag = 3

		(while the first algorithm is app. 3 times faster than the third one, 
		the third algorithm should be in general much more precise).  
		The computed frequencies are in the range given by minfreq and maxfreq.
		The function returns the number of determined frequencies or 0 in the case
		of error.

		The vectors input[1][j] and input[2][j], j = 1 ... ndata (ndata must
		be a power of 2), are the input data X(j-1) and Y(j-1).
		"""
		output = np.empty((nfreq,3),order='C',dtype=np.float64)
		try:
			self._fmft( output, nfreq, minfreq, maxfreq,flag, np.array(inpt,order='C',dtype=np.float64)  , c_long(ndata) )
			return output
		except RuntimeError:
			print "FMFT failed!"
			return output

	def fmft(self,nfreq,inpt):
                """
                A simplified version of the full fmft routine that automatically computes
                the number of data points to use based on the length of the input data and
                limits the frequency range between -1/2 to 1/2 cycles per sampling period,
                i.e., the Nyquist frequencies.
                """
                ndata = nearest_pow2(len(inpt))
                nyq = 2*np.pi * 0.49
                method = 3
                return self.fmft_full(nfreq,-1.*nyq,nyq,method,inpt,ndata)

def fmft_plot(wrap,nfreq,inpt):
	
	out = wrap.fmft(nfreq,inpt)
	times = inpt[:,0]
	# Get FMFT-approximation of data
	simdata  = np.zeros((len(times),2))
	for o in out:
		w,a,p = o
		simdata[:,0] = simdata[:,0] + a*np.cos(w * times + p)
		simdata[:,1] = simdata[:,1] + a*np.sin(w * times + p) 
	plt.plot(times,simdata[:,0],'b--')
	plt.plot(times,simdata[:,1],'g--')			
	plt.plot(times,inpt[:,1],'k-')
	plt.plot(times,inpt[:,2],'r-')			
	plt.show()
if __name__=="__main__":

	npts = 2*512
	times = np.linspace(0,2*np.pi*300,npts)
	data = np.zeros((npts,2))

	for w in [0.1 ,0.01, 0.03, 0.02]:
		
		amp = np.random.rand()
		phase = np.random.rand()
		data[:,0] = data[:,0] + amp*np.cos(w*times + phase) 
		data[:,1] = data[:,1] + amp*np.sin(w*times + phase)
		
	figure()
	plot(times,data[:,0],'k-')
	plot(times,data[:,1],'r-')			
	
	ff = libwrapper()
	fulldata = np.hstack((times.reshape(-1,1),data))
	out = ff.fmft(3,fulldata)

	simdata  = np.zeros((npts,2))
	for o in out:
		print o
		w,a,p = o
		simdata[:,0] = simdata[:,0] + a*np.cos(w * times + p)
		simdata[:,1] = simdata[:,1] + a*np.sin(w * times + p) 
		
	plot(times,simdata[:,0],'y--')
	plot(times,simdata[:,1],'g--')			
	show()

