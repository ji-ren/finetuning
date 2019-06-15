import time
start=time.clock()
import sympy
import FT_jxl_analytical

#set up parameters
mh2=sympy.Symbol('mh2')		
g1=sympy.Symbol('g1')		
g2=sympy.Symbol('g2')		
G=g1**2+g2**2
mt1=sympy.Symbol('mt1')		
mt2=sympy.Symbol('mt2')		
mt=sympy.Symbol('mt')		
At=sympy.Symbol('At')		
MHU2=sympy.Symbol('MHU2')	
MHD2=sympy.Symbol('MHD2')	
MS2=sympy.Symbol('MS2')		
A=sympy.log(mt1*mt2/mt**2)
B=At/(mt1**2-mt2**2)*sympy.log(mt1**2/mt2**2)
C=At**2/(mt1**2-mt2**2)**2*(1-(mt1**2+mt2**2)/(mt1**2-mt2**2)*sympy.log(mt1/mt2))
TB=sympy.Symbol('TB')
TB2=TB**2
Mz=sympy.Symbol('Mz')
MU=sympy.Symbol('MU')
AL=sympy.Symbol('AL')
AK=sympy.Symbol('AK')
L=sympy.Symbol('L')	
K=sympy.Symbol('K')	
yt=sympy.Symbol('yt')
origin=3*yt**4/(8*(sympy.pi)**2)*2*Mz**2*TB2/(G*(1+TB2))
Wt=origin*(A+B*(At+MU/TB)+C*(At+MU/TB)**2)

#获得数值
V_mh2=15703.778547419524
V_g1=float(3.61496259E-01)   #spc.Value('GAUGE',[1])
V_g2=float(6.71344778E-01)   #spc.Value('GAUGE',[2])
V_mt1=float(2.01037642E+03)   #spc.Value('MASS',[1000006])
V_mt2=float(2.12820584E+03)   #spc.Value('MASS',[2000006])
V_mt=float(1.73500000E+02)    #spc.Value('SMINPUTS',[6])
V_At=2000.0    #spc.Value('EXTPAR',[11])
V_MHU2=float(-3.40835683E+05)     #spc.Value('MSOFT',[22])
V_MHD2=float(3.45021908E+06)     #spc.Value('MSOFT',[21])
V_MS2=float(1.79817995E+04)      #spc.Value('NMSSMRUN',[10])
V_TB=float(1.14531792E+01)    #spc.Value('MINPAR',[3])
V_Mz=float(9.11887000E+01)    #spc.Value('MASS',[23])
V_MU=300.0     #spc.Value('NMSSMRUN',[5])
V_AL=2000.0    #spc.Value('NMSSMRUN',[3])
V_AK=200.0   #spc.Value('NMSSMRUN',[4])
V_L=float(2.28888827E-01)    #spc.Value('NMSSMRUN',[1])
V_K=float(-1.39604166E-01)   #spc.Value('NMSSMRUN',[2])
V_yt=0.003        #spc.Value('YU',[3])
#V_yt=2**0.5*V_mt/(sympy.sin(sympy.atan(V_TB))*246)
#V_yt=V_mt/(sympy.sin(sympy.atan(V_TB))*174)

V_dic={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt}
fine_tuning=FT_jxl_analytical.fine_tuning
#fine_tuning=finetuning(Mz_FT1,Mz_FT2,Mz_FT3,Mz_FT4,Mz_FT5,Mz_FT6,Mz_FT7,Mz_FT8,mh_FT1,mh_FT2,mh_FT3,mh_FT4,mh_FT5,mh_FT6,mh_FT7,mh_FT8)
print(fine_tuning.Mz_FT1.evalf(subs=V_dic))
print(fine_tuning.Mz_FT2.evalf(subs=V_dic))
print(fine_tuning.Mz_FT3.evalf(subs=V_dic))
print(fine_tuning.Mz_FT4.evalf(subs=V_dic))
print(fine_tuning.Mz_FT5.evalf(subs=V_dic))
print(fine_tuning.Mz_FT6.evalf(subs=V_dic))
print(fine_tuning.Mz_FT7.evalf(subs=V_dic))
print(fine_tuning.Mz_FT8.evalf(subs=V_dic))


end=time.clock()
print('time:',end-start)
