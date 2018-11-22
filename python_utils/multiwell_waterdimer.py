import numpy as np

beta1=12000
beta2=20000

def calcdelta(Ia,Ib):
    y1= np.arctanh(Ia)
    y2= np.arctanh(Ib)
    delta= 2.0*(y2-y1)/(beta2-beta1)
    return(delta)

def calcbetabar(Ia, Ib):
    y1= np.arctanh(Ia)
    y2= np.arctanh(Ib)
    betabar= (beta1*y2 - beta2*y1)/(y2-y1)
    return(betabar)

n=8

c= np.zeros((6,8))
#E, A, D, G, G, AG, AG, B
#E, (12), (34), (ab)(13)(24), (ab)(14)(23), (ab)(1324), (ab)(1423), (12)(34),

c[0,:]= (1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)#/np.sqrt(8.0) #A1+
c[1,:]= (1.0,1.0,-1.0,0.0,0.0,0.0,0.0,-1.0)#/np.sqrt(4.0) #E+
c[2,:]= (1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0)#/np.sqrt(8.0) #B1+
c[3,:]= (1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0)#/np.sqrt(8.0) #A2-
c[4,:]= (1.0,-1.0,1.0,0.0,0.0,0.0,0.0,-1.0)#/np.sqrt(4.0) #E-
c[5,:]= (1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0)#/np.sqrt(8.0)#B2-
# print np.shape(c)
# print c

A1=2.0*0.111478983128910#0.112884788776427
A2=2.0*0.230546268494811#0.230546268494811
AG1=1.071538742584852E-003#7.699590610213557E-004
AG2=3.118506023220292E-003#3.312585739684554E-003
G1=5.632020635958411E-003
G2=1.049548502801190E-002
B1=0.0#1.056772155837994E-006*4.0
B2=0.0#4.846535050718374E-006*4.0
D1=0.0#1.900306551351940E-005
D2=0.0#3.773895955538238E-004

# A1=2.0*1.340366910821415E-002
# A2=2.0*0.102086996302219
# AG1=7.699590610213557E-004
# AG2=3.312585739684554E-003
# G1=4.761792525321749E-003
# G2=1.108151783073477E-002
# B1=1.056772155837994E-006*4.0
# B2=4.846535050718374E-006*4.0
# D1=1.900306551351940E-005
# D2=3.773895955538238E-004
#E, A, D, AG, AG, G, G, B
#E, (12), (34), (ab)(13)(24), (ab)(14)(23), (ab)(1324), (ab)(1423), (12)(34),

rho1=np.asarray([1.0, A1, D1, AG1, AG1, G1, G1, B1])
rho2=np.asarray([1.0, A2, D2, AG2, AG2, G2, G2, B2])

etaminus1= np.zeros((n))
etaplus1= np.zeros((n))
etaminus2= np.zeros((n))
etaplus2= np.zeros((n))

I1= np.zeros(6)
I2= np.zeros(6)
for i in range(6):
    for j in range(n):
        etaminus1[i]+= (1.0- c[i,j])*rho1[j]
        etaplus1[i]+= (1.0+ c[i,j])*rho1[j]
        etaminus2[i]+= (1.0- c[i,j])*rho2[j]
        etaplus2[i]+= (1.0+ c[i,j])*rho2[j]
    I1[i]= etaminus1[i]/etaplus1[i]
    I2[i]= etaminus2[i]/etaplus2[i]

labels= ['A1+', 'E+', 'B1+', 'A2-', 'E-', 'B2-']
levels= np.zeros_like(etaminus1)
betabars=np.zeros_like(etaminus1)
outputfile= open("pimd.dat","w")
for i in range(6):
    levels[i]=219475.0*calcdelta(I1[i],I2[i])
    betabars[i]= calcbetabar(I1[i], I2[i])
    print labels[i],"=", levels[i], betabars[i], I1[i],I2[i]
    outputfile.write(str(levels[i]) + "\n")

# h17= (levels[3] - levels[4] + levels[1])/6.0
# h15= (levels[1] - 4.0*h17)/4.0
# h13= (levels[4] - levels[5] + levels[2] - 4.0*h15)/16.0
# h14= 2.0*h13 + (levels[5] - levels[2])/4.0
# h12= 0.5*(levels[2] - 8.0*h13 - 2.0*h15 - 2.0*h17)

# print "------------------"
# print "h12=", h12
# print "h13=", h13
# print "h14=", h14
# print "h15=", h15
# print "h17=", h17
# print "------------------"
# print "Splittings:"
# print "acceptor=", 4.0*h14
# print "lower interchange=", 4.0*(h15 + h17)
# print "upper interchange=", 4.0*(h15-h17)
# print "lower bifurcation=", h12 + 4.0*h13
# print "upper bifurcation=", h12 - 4.0*h13
