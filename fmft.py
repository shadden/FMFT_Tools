from ctypes import *
import numpy as np

LIBPATH = "/Users/samuelhadden/19_FMFT"

def get_ctype_ptr(dtype,dim,**kwargs):
	return np.ctypeslib.ndpointer(dtype=dtype,ndim=dim,flags='C',**kwargs)
	
p1d=get_ctype_ptr(np.float,1)
p2d = get_ctype_ptr(np.float,2)
pd = POINTER(c_double)
ppd = POINTER(pd)
p1dInt = get_ctype_ptr(np.int,1)

def nearest_pow2(x):
	int(2**np.floor(np.log2(x)))
	
class libwrapper(object):
	""" A wrapper class for the FMFT library.  Gives access to the function fmft"""
	def __init__(self):
		
		self.lib = CDLL("%s/libfmft.so" % LIBPATH)
		self._fmft = self.lib.fmftWrap
		self._fmft.argtypes =[p2d, c_int,c_double, c_double, c_int, p2d, c_long]
		self._fmft.restype = c_int
		def check_errors(ret, func, args):
			if ret<=0:
				raise RuntimeError("FMFT returned error code %d for the given arguments"%ret)
			return ret
		self._fmft.errcheck = check_errors
	def fmft(self, nfreq, minfreq, maxfreq, flag, input, ndata):
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
		output = np.empty((nfreq,3),order='C')
		try:
			self._fmft( output, nfreq, minfreq, maxfreq,flag, np.array(input,order='C') , ndata)
			return output
		except RuntimeError:
			print "FMFT failed!"
			return output
if __name__=="__main__":
	ff = libwrapper()
	x = np.loadtxt("tmp/tp177.aei",skiprows=4)
	zz=np.vstack(( x[:,0], x[:,2]*np.cos(x[:,4]*np.pi/180.),x[:,2]*np.sin(x[:,4]*np.pi/180.) )).T
	out = ff.fmft(4,-.3,.3,2,zz,128)
	for l in out:
		print l