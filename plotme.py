from argparse import ArgumentParser

parser = ArgumentParser(description='plot N-body output')
parser.add_argument('-d','--dir',metavar='DIR',default="Integration",help='Source directory of *.aei files')
args = parser.parse_args()
dir = args.dir

def slope(t,val):
	interp,s = linalg.lstsq( vstack([ones(len(t)),t]).T , val )[0]
	return s

t,a,e,i,per,node,M,mass,x,y,z,vx,vy,vz=loadtxt("%s/PL1.aei"%dir,skiprows=4).T
t,a1,e1,i1,per1,node1,M1,mass1,x1,y1,z1,vx1,vy1,vz1=loadtxt("%s/PL2.aei"%dir,skiprows=4).T


pmg =  (per + node) * pi /180.
pmg1 = (per1 + node1) * pi / 180.

ex,ey = e*cos(pmg), e*sin(pmg)
ex1,ey1 = e1*cos(pmg1), e1*sin(pmg1)

l1 =(per1 + node1 + M1) * pi / 180. 
l = (per + node + M) * pi / 180.
phi = 3 * l1 - 2 * l - pmg
phi1= 3 * l1 - 2 * l - pmg1

phi2 = 6 * l1 - 4 * l - pmg - pmg1

phi,phi1 = mod(phi,2*pi),mod(phi1,2*pi)
phi2 = mod(phi2,2*pi)



l0,n = linalg.lstsq( vstack([ones(len(t)),t]).T , unwrap(l) )[0]
l10,n1 = linalg.lstsq( vstack([ones(len(t)),t]).T , unwrap(l1) )[0]
dl = unwrap(l) - l0 - n*t 
dl1 = unwrap(l1) - l10 - n1*t 

figure()
##
ax=subplot(221)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 
elim=max(map(max,[abs(ey),abs(ex)]))
plot(ex,ey)
xlim(-1*elim,elim)
ylim(-1*elim,elim)
ax.set_xticks(ax.get_xticks()[::2])
ax.set_yticks(ax.get_yticks()[1::2])

##
ax=subplot(222)
ax.get_yaxis().set_major_formatter(NullFormatter())
plot(t,ey,'r')
ylim(-1*elim,elim)
ax.set_xticks(ax.get_xticks()[1:])
##
ax=subplot(223)
plot(ex[::-1],t[::-1],'g')
ax.set_yticks(ax.get_yticks()[1:-1])
ax.get_xaxis().set_major_formatter(NullFormatter())
xlim(-1*elim,elim)
#####
subplots_adjust(hspace=0.0,wspace=0.0)
####
########
figure()
x,y = e*sin(-pmg), e*cos(-pmg)
zx = y * cos(3*l1 - 2*l) - x * sin(3*l1 - 2*l)
zy = x * cos(3*l1 - 2*l) + y * sin(3*l1 - 2*l)
plot(zx,zy,'r-')
x1,y1 = e1*sin(-pmg1), e1*cos(-pmg1)
zx = y1 * cos(3*l1 - 2*l) - x1 * sin(3*l1 - 2*l)
zy = x1 * cos(3*l1 - 2*l) + y1 * sin(3*l1 - 2*l)
plot(zx,zy,'b-')

title('e * exp{i(3l1 - 2l - w)}')

figure()
subplot(211)
fmftdat = loadtxt("fmftA.out",skiprows=1);
amp1 = fmftdat[1,1]
phi1 = fmftdat[1,2] * pi /180.
freq1 = fmftdat[1,0]
plot(t, 2* amp1* cos( freq1 *t + phi1) - 2* amp1* cos( freq1 * 0 + phi1) ,'k.')
plot(t,a/a[0] - 1,'b-')
xlim(xmax=500)

subplot(212)
###
fmftdat = loadtxt("fmftL.out",skiprows=1);
amp1 = fmftdat[0,1]
phi1 = fmftdat[0,2] * pi /180.
freq1 = fmftdat[0,0]
##
amp2 = fmftdat[2,1]
phi2 = fmftdat[2,2] * pi /180.
freq2 = fmftdat[2,0]
###
plot(t, 2* amp1* cos( freq1 *t + phi1) + 2* amp2* cos( freq2 *t + phi2) ,'k.')
plot(t,dl,'r-')
xlim(xmax=500)

######################################################
alpha = -2.02522 
beta = 2.48401
C=sqrt(alpha**2 + beta**2)
v1 = (alpha * ey + beta * ey1)/C
u1 = (-1*(alpha * ex + beta * ex1))/C
u2 = (beta * (-1.*ex) - alpha * (-1.*ex1))/C
v2 = (beta * ey - alpha * ey1)/C
P1 = 0.5 * (u1**2 + v1**2)
P2 = 0.5 * (u2**2 + v2**2)

psi1 = mod(3*l1 - 2 * l + arctan(v1,u1),2*pi)
psi2 = mod(3*l1 - 2 * l + arctan(v2,u2),2*pi)
figure()
tlim=1.e3
plot(sqrt(2*P1[t<tlim])*cos(psi1[t<tlim]),sqrt(2*P1[t<tlim])*sin(psi1[t<tlim]),'k-')
title('Z * exp{i(3l1 - 2l)}')
######################################################
show()
