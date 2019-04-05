import numpy as np
import math
from numpy.linalg import inv
from numpy import linalg as LA
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm
ath=10.456788 
ath=float(ath)
btha=10.456788 
btha=float(btha)
cth=11.769562 
c=float(cth)
n=0
n=float(n)
alphath=90.000000 
alphath=float(alphath)
betath=90.000000
betath=float(betath)
gammath=90.000000
gammath=float(gammath)
aex=10.56720
aec=float(aex)
bexa=10.56720
bexa=float(bexa)
cex=11.9841
cex=float(cex)
alphaex=90.00
alphaex=float(alphaex)
betaex=90.000000 
betaex=float(betaex)
gammaex=90.00 
gammaex=float(gammaex)
ascal=1.0

 
y = sum(1 for line in open('inpex.txt'))
for x in range (0,y):

 arraysex = [np.array(map(str,line.split())) for line in open('inpex.txt')]
 arraysth = [np.array(map(str,line.split())) for line in open('inpth.txt')]
 arraysex=np.append(arraysex[x],[0,0,0])
 arraysth=np.append(arraysth[x],[0,0,0])
 bex = np.reshape(arraysex, (2, 5))
 bth = np.reshape(arraysth, (2, 5))
 eex=np.delete(bex,[0,0])
 eth=np.delete(bth,[0,0])
 newex=np.reshape(eex,[3,3])
 newth=np.reshape(eth,[3,3])
 Aex=np.array([[float(newex[0,0]), float(newex[1,2]),float(newex[1,1])], [float(newex[1,2]),float(newex[0,1]),float(newex[1,0])],[float(newex[1,1]),float(newex[1,0]),float(newex[0,2])]])
 Ath=np.array([[ascal*float(newth[0,0]), ascal*float(newth[1,2]),ascal*float(newth[1,1])], [ascal*float(newth[1,2]),ascal*float(newth[0,1]),ascal*float(newth[1,0])],[ascal*float(newth[1,1]),ascal*float(newth[1,0]),ascal*float(newth[0,2])]])
 B2th=np.array([[ath,btha*(np.cos(gammath*np.pi/180.0)),cth*np.cos(betath*np.pi/180.0)],[n,btha*np.sin(gammath*np.pi/180.0),cth*((np.cos(alphath*np.pi/180.0)-(np.cos(betath*np.pi/180.0))*np.cos(gammath*np.pi/180.0)))/np.sin(gammath*np.pi/180.0)],[n,n,(ath*btha*cth*np.sqrt(1-np.cos(alphath*np.pi/180.0)*np.cos(alphath*np.pi/180.0)-np.cos(betath*np.pi/180.0)*np.cos(betath*np.pi/180.0)-np.cos(gammath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)+2*np.cos(alphath*np.pi/180.0)*np.cos(betath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)))/(ath*btha*np.sin(gammath*np.pi/180.0))]])
 B2ex=np.array([[aex,bexa*(np.cos(gammaex*np.pi/180.0)),cex*np.cos(betaex*np.pi/180.0)],[n,bexa*np.sin(gammaex*np.pi/180.0),cex*((np.cos(alphaex*np.pi/180.0)-(np.cos(betaex*np.pi/180.0))*np.cos(gammaex*np.pi/180.0)))/np.sin(gammaex*np.pi/180.0)],[n,n,(aex*bexa*cex*np.sqrt(1-np.cos(alphaex*np.pi/180.0)*np.cos(alphaex*np.pi/180.0)-np.cos(betaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)-np.cos(gammaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)+2*np.cos(alphaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)))/(aex*bexa*np.sin(gammaex*np.pi/180.0))]])
 Mth=np.transpose(B2th)
 Mex=np.transpose(B2ex)
 Aex2=np.matmul(B2ex,Aex)
 Aex22=np.matmul(Aex2,Mex)
 Ath2=np.matmul(B2th,Ath)
 Ath22=np.matmul(Ath2,Mth)
 wexpv = LA.eigvals(Aex22)
 wexps=np.sort(wexpv)
 wthpv = LA.eigvals(Ath22)
 wthps=np.sort(wthpv)
 f = open('outexp.txt', 'a')
 print>>f,bex[0,0],wexps[0],wexps[1],wexps[2]
 g = open('outth.txt','a')
 print>>g,bth[0,0],wthps[0],wthps[1],wthps[2]
 Aexmin=np.amin([[wexpv[0],wexpv[1],wexpv[2]]])
 Aexmax=np.amax([[wexpv[0],wexpv[1],wexpv[2]]])
 Aw=Aexmax/Aexmin

 Aex=inv(Aex22)
 Ath=inv(Ath22) 
 wex,vex = LA.eig(Aex)
 wth,vth = LA.eig(Ath)
 uex1 = np.array([vex[0,0], vex[1,0], vex[2,0] ])
 uex2 = np.array([vex[0,1], vex[1,1], vex[2,1] ])
 uex3 = np.array([vex[0,2], vex[1,2], vex[2,2] ])
 uth1 = np.array([vth[0,0], vth[1,0], vth[2,0] ])
 uth2 = np.array([vth[0,1], vth[1,1], vth[2,1] ])
 uth3 = np.array([vth[0,2], vth[1,2], vth[2,2] ])
 c1 = dot(uex1,uth1)/norm(uex1)/norm(uth1) 
 c2 = dot(uex1,uth2)/norm(uex1)/norm(uth2)
 c3 = dot(uex1,uth3)/norm(uex1)/norm(uth3)
 c4 = dot(uex2,uth1)/norm(uex2)/norm(uth1)
 c5 = dot(uex2,uth2)/norm(uex2)/norm(uth2)
 c6 = dot(uex2,uth3)/norm(uex2)/norm(uth3)
 c7 = dot(uex3,uth1)/norm(uex3)/norm(uth1)
 c8 = dot(uex3,uth2)/norm(uex3)/norm(uth2)
 c9 = dot(uex3,uth3)/norm(uex3)/norm(uth3)
 angle1 = arccos(clip(c1, -1, 1)) 
 angle2 = arccos(clip(c2, -1, 1))
 angle3 = arccos(clip(c3, -1, 1))
 angle4 = arccos(clip(c4, -1, 1)) 
 angle5 = arccos(clip(c5, -1, 1)) 
 angle6 = arccos(clip(c6, -1, 1)) 
 angle7 = arccos(clip(c7, -1, 1)) 
 angle8 = arccos(clip(c8, -1, 1)) 
 angle9 = arccos(clip(c9, -1, 1)) 
 Volth=(1/(math.sqrt(abs(wth[0]))))*(1/(math.sqrt(abs(wth[1]))))*(1/(math.sqrt(abs(wth[2]))))
 Volex=(1/(math.sqrt(abs(wex[0]))))*(1/(math.sqrt(abs(wex[1]))))*(1/(math.sqrt(abs(wex[2]))))
 Volver=Volth/Volex
 Vectorth=np.array([[1/ath,-(np.cos(gammath*np.pi/180.0))/(ath*np.sin(gammath*np.pi/180.0)),(btha*cth*(np.cos(alphath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)-np.cos(betath*np.pi/180.0)))/(np.sin(gammath*np.pi/180.0)*ath*btha*cth*np.sqrt(1-np.cos(alphath*np.pi/180.0)*np.cos(alphath*np.pi/180.0)-np.cos(betath*np.pi/180.0)*np.cos(betath*np.pi/180.0)-np.cos(gammath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)+2*np.cos(alphath*np.pi/180.0)*np.cos(betath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)))],[0,1/(btha*np.sin(gammath*np.pi/180.0)),ath*cth*(np.cos(betath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)-np.cos(alphath*np.pi/180.0))/(np.sin(gammath*np.pi/180.0)*ath*btha*cth*np.sqrt(1-np.cos(alphath*np.pi/180.0)*np.cos(alphath*np.pi/180.0)-np.cos(betath*np.pi/180.0)*np.cos(betath*np.pi/180.0)-np.cos(gammath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)+2*np.cos(alphath*np.pi/180.0)*np.cos(betath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)))],[0,0,(ath*btha*np.sin(gammath*np.pi/180.0))/(ath*btha*cth*np.sqrt(1-np.cos(alphath*np.pi/180.0)*np.cos(alphath*np.pi/180.0)-np.cos(betath*np.pi/180.0)*np.cos(betath*np.pi/180.0)-np.cos(gammath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)+2*np.cos(alphath*np.pi/180.0)*np.cos(betath*np.pi/180.0)*np.cos(gammath*np.pi/180.0)))]])
 Vectorexp=np.array([[1/aex,-(np.cos(gammaex*np.pi/180.0))/(aex*np.sin(gammaex*np.pi/180.0)),(bexa*cth*(np.cos(alphaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)-np.cos(betaex*np.pi/180.0)))/(np.sin(gammaex*np.pi/180.0)*aex*bexa*cex*np.sqrt(1-np.cos(alphaex*np.pi/180.0)*np.cos(alphaex*np.pi/180.0)-np.cos(betaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)-np.cos(gammaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)+2*np.cos(alphaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)))],[0,1/(bexa*np.sin(gammaex*np.pi/180.0)),aex*cex*(np.cos(betaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)-np.cos(alphaex*np.pi/180.0))/(np.sin(gammaex*np.pi/180.0)*aex*bexa*cex*np.sqrt(1-np.cos(alphaex*np.pi/180.0)*np.cos(alphaex*np.pi/180.0)-np.cos(betaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)-np.cos(gammaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)+2*np.cos(alphaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)))],[0,0,(aex*bexa*np.sin(gammaex*np.pi/180.0))/(aex*bexa*cex*np.sqrt(1-np.cos(alphaex*np.pi/180.0)*np.cos(alphaex*np.pi/180.0)-np.cos(betaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)-np.cos(gammaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)+2*np.cos(alphaex*np.pi/180.0)*np.cos(betaex*np.pi/180.0)*np.cos(gammaex*np.pi/180.0)))]])
 h = open('out.txt', 'a')
 if Aexmax/Aexmin>1.6:
   print>>h, "Anisotrop", bex[0,0], Aexmax/Aexmin,np.matmul(Vectorexp,uex1),wex[0], (np.matmul(Vectorexp,uex2)),wex[1],(np.matmul(Vectorexp,uex3)),wex[2],(np.matmul(Vectorth,uth1)),wth[0],(np.matmul(Vectorth,uth2)),wth[1],(np.matmul(Vectorth,uth3)),wth[2]
 print>>h, bex[0,0],"Verhaltnis",Volver,"Winkel",  "(11)",angle1*57.296, "(22)",angle5*57.296, "(33)",angle9*57.296, "(12)",angle2*57.296, "(21)",angle4*57.296, "(33)",angle9*57.296, "(13)",angle3*57.296, "(21)",angle4*57.296%90, "(32)",angle8*57.296
 print>>h, "(11)",  angle1*57.296, "(23)", angle6*57.296, "(32)", angle8*57.296, "(12)", angle2*57.296, "(23)",angle6*57.296, "(31)", angle7*57.296, "(13)", angle3*57.296, "(22)",angle5*57.296, "(31)", angle7*57.296
 
 
 


 

 
 
 
 


    
    
     
