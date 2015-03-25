from numpy import *
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description='run FMFT on N-body output')
parser.add_argument('-s','--semimajorA',default=False,action='store_true',help='Include analysis of da/a')
parser.add_argument('--dir',metavar='DIR',default="Integration",help='Source directory of *.aei files')
parser.add_argument('-t','--testparticle',default=False,action='store_true',help='Test particle: ignore outer planet files')

args = parser.parse_args()
tp = args.testparticle
#########################################################################
dir = args.dir
data = loadtxt("%s/PL1.aei"%dir,skiprows=4)
data1 = loadtxt("%s/PL2.aei"%dir,skiprows=4)
t = data.T[0]
t1 = data1.T[0]
dt = t[1] - t[0]
dt1 = t1[1] - t1[0]

l = unwrap( (data.T[6] + data.T[4] + data.T[5]) * pi / 180. )
l1 = unwrap( ( data1.T[6] + data1.T[4] + data1.T[5] ) * pi / 180. )

x,y = data.T[8],data.T[9]
x1,y1 = data1.T[8],data1.T[9]
theta = unwrap( arctan2(y,x) )
theta1 = unwrap( arctan2(y1,x1) )

l0,n = linalg.lstsq( vstack([ones(len(t)),t]).T , l )[0]
l10,n1 = linalg.lstsq( vstack([ones(len(t)),t]).T , l1 )[0]

dl = l - l0 - n*t 
dl1 = l1 - l10 - n1*t 

theta0,ptheta = linalg.lstsq( vstack([ones(len(t)),t]).T , theta )[0]
theta10,ptheta1 = linalg.lstsq( vstack([ones(len(t)),t]).T , theta1 )[0]

dtheta = theta - theta0 - ptheta * t
dtheta1 = theta1 - theta10 - ptheta1 * t

#########################################################################
with open("outputE.dat","w") as fi:
	for d in data:
		ecc = d[2]
		pmg = (d[4]+d[5]) * pi/180.
		fi.write( "%.12f\t%.12f\t%.12f\n" % (d[0], ecc * cos(pmg), ecc * sin(pmg)) )
with open("outputE1.dat","w") as fi:
	for d in data1:
		ecc = d[2]
		pmg = (d[4]+d[5]) * pi/180.
		fi.write( "%.12f\t%.12f\t%.12f\n" % (d[0], ecc * cos(pmg), ecc * sin(pmg)) )

commandstring = "./samfmft %f < %s" % (dt,"outputE.dat")

proc = subprocess.Popen(commandstring,\
        stdout=subprocess.PIPE,shell=True)

(out,err) = proc.communicate()

print "Inner planet eccentricity: "
print out
with open("fmftE.out","w") as fi:
	fi.writelines(out)

if not tp:
	commandstring = "./samfmft %f < %s" % (dt1,"outputE1.dat")

	proc = subprocess.Popen(commandstring,\
			stdout=subprocess.PIPE,shell=True)

	(out,err) = proc.communicate()

	print "Outer planet eccentricity: "
	print out
	with open("fmftE1.out","w") as fi:
		fi.writelines(out)

#########################################################################
with open("outputL.dat","w") as fi:
	for d in zip(t,dl):
		fi.write( "%.12f\t%.12f\t%.12f\n" % (d[0], d[1] , 0.0) )
with open("outputL1.dat","w") as fi:
	for d in zip(t1,dl1):
		fi.write( "%.12f\t%.12f\t%.12f\n" % (d[0], d[1] , 0.0) )

commandstring = "./samfmft %f < %s" % (dt,"outputL.dat")

proc = subprocess.Popen(commandstring,\
        stdout=subprocess.PIPE,shell=True)

(out,err) = proc.communicate()

print "Inner planet lambda: "
print out
with open("fmftL.out","w") as fi:
	fi.writelines(out)

if not tp:
	commandstring = "./samfmft %f < %s" % (dt1,"outputL1.dat")

	proc = subprocess.Popen(commandstring,\
			stdout=subprocess.PIPE,shell=True)

	(out,err) = proc.communicate()

	print "Outer planet lambda: "
	print out
	with open("fmftL1.out","w") as fi:
		fi.writelines(out)

#########################################################################
with open("outputTH.dat","w") as fi:
	for d in zip(t,dtheta):
		fi.write( "%.5f\t%.5f\t%.5f\n" % (d[0], d[1] , 0.0) )
with open("outputTH1.dat","w") as fi:
	for d in zip(t1,dtheta1):
		fi.write( "%.5f\t%.5f\t%.5f\n" % (d[0], d[1] , 0.0) )

commandstring = "./samfmft %f < %s" % (dt,"outputTH.dat")

proc = subprocess.Popen(commandstring,\
        stdout=subprocess.PIPE,shell=True)

(out,err) = proc.communicate()

print "Inner planet theta: "
print out
with open("fmftTH.out","w") as fi:
	fi.writelines(out)

if not tp:
	commandstring = "./samfmft %f < %s" % (dt1,"outputTH1.dat")

	proc = subprocess.Popen(commandstring,\
			stdout=subprocess.PIPE,shell=True)

	(out,err) = proc.communicate()

	print "Outer planet theta: "
	print out
	with open("fmftTH1.out","w") as fi:
		fi.writelines(out)

#########################################################################
if args.semimajorA:

	daBya = (data.T[1] / data[0,1]) -1.
	da1Bya1 = (data1.T[1] / data1[0,1]) -1.

	with open("outputA.dat","w") as fi:
		for d in zip(t,daBya):
			fi.write( "%.8f\t%.8f\t%.8f\n" % (d[0], d[1] , 0.0) )
	with open("outputA1.dat","w") as fi:
		for d in zip(t1,da1Bya1):
			fi.write( "%.5f\t%.5f\t%.5f\n" % (d[0], d[1] , 0.0) )
	
	commandstring = "./samfmft %f < %s" % (dt,"outputA.dat")
	
	proc = subprocess.Popen(commandstring,\
	        stdout=subprocess.PIPE,shell=True)
	
	(out,err) = proc.communicate()
	
	print "Inner planet semi-major axis: "
	print out
	with open("fmftA.out","w") as fi:
		fi.writelines(out)

	if not tp:
		commandstring = "./samfmft %f < %s" % (dt1,"outputA1.dat")

		proc = subprocess.Popen(commandstring,\
				stdout=subprocess.PIPE,shell=True)

		(out,err) = proc.communicate()

		commandstring = "./samfmft %f < %s" % (dt1,"outputA1.dat")
	
		proc = subprocess.Popen(commandstring,\
				stdout=subprocess.PIPE,shell=True)
	
		(out,err) = proc.communicate()
	
		print "Outer planet semi-major axis: "
		print out
		with open("fmftA1.out","w") as fi:
			fi.writelines(out)

