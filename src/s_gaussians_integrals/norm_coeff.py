#Module for computing normalization coefficients for an arbitrary
#uncontracted s-Gaussian

def n_cof(alpha):
	import math
	#print math.pi
	return (math.pi/(2.0*alpha))**(-3.0/4.0)

def overlap(alp1,x1,y1,z1,alp2,x2,y2,z2):
	import math
	r12 = math.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2 )
	return n_cof(alp1)*n_cof(alp2)*math.exp(-alp1*alp2/(alp1+alp2)*r12**2)\
	*(math.pi/(alp1+alp2))**(3.0/2.0)

def kin(alp1,x1,y1,z1,alp2,x2,y2,z2):
	import math
	r12 = math.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2 )
	return n_cof(alp1)*n_cof(alp2)*(-2.0*(alp1*alp2/(alp1+alp2))**2*r12**2\
	+ 3*(1.0*alp1*alp2/(alp1+alp2)))*math.exp(-1.0*alp1*alp2/(alp1+alp2)*r12**2)*\
	(math.pi/(alp1+alp2))**(3.0/2.0)

def nuc(alp1,x1,y1,z1,alp2,x2,y2,z2,Zn,nx,ny,nz):
	import math 
	r12 = math.sqrt((float(x1)-float(x2))**2+(float(y1)-float(y2))**2+(float(z1)-float(z2))**2 )
	xc = (1.0*alp1*x1+1.0*alp2*x2)/(alp1+alp2)
	yc = (1.0*alp1*y1+1.0*alp2*y2)/(alp1+alp2)
	zc = (1.0*alp1*z1+1.0*alp2*z2)/(alp1+alp2)

	rac = math.sqrt((float(xc)-float(nx))**2+(float(yc)-float(ny))**2+(float(zc)-float(nz))**2 ) 
	if rac == 0.0:
		return -2.0*Zn*n_cof(alp1)*n_cof(alp2)*(math.pi/(alp1+alp2))**(3.0/2.0)*\
			math.exp(-1.0*alp1*alp2/(alp1+alp2)*r12**2)
	elif rac !=0.0: 
		return -1.0*Zn*n_cof(alp1)*n_cof(alp2)*(math.pi/(alp1+alp2))**(3.0/2.0)*\
                        math.exp(-1.0*alp1*alp2/(alp1+alp2)*r12**2)*\
			math.erf(math.sqrt(alp1+alp2)*rac)/rac

