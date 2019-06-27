import time
start=time.clock()
import xslha
import sympy 

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
Mz2=sympy.Symbol('Mz2')
MU=sympy.Symbol('MU')
AL=sympy.Symbol('AL')
AK=sympy.Symbol('AK')
L=sympy.Symbol('L')	
K=sympy.Symbol('K')	
yt=sympy.Symbol('yt')
origin=3*yt**4/(8*(sympy.pi)**2)*2*Mz2*TB2/(G*(1+TB2))
#Wt=origin*sympy.log(mt1**2/mt2**2)
Wt=origin*(A+B*(At+MU/TB)+C*(At+MU/TB)**2)
# print('导入参数完成')

# #获得数值
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
V_Mz2=float(9.11887000E+01)**2    #spc.Value('MASS',[23])
V_MU=300.0     #spc.Value('NMSSMRUN',[5])
V_AL=2000.0    #spc.Value('NMSSMRUN',[3])
V_AK=200.0   #spc.Value('NMSSMRUN',[4])
V_L=float(2.28888827E-01)    #spc.Value('NMSSMRUN',[1])
V_K=float(-1.39604166E-01)   #spc.Value('NMSSMRUN',[2])
V_yt=0       #spc.Value('YU',[3])
#V_yt=2**0.5*V_mt/(sympy.sin(sympy.atan(V_TB))*246)
#V_yt=V_mt/(sympy.sin(sympy.atan(V_TB))*174)
# print('获取数值完成')

# spc=xslha.read('/home/jxl/Desktop/tools/NMSSMTools_5.4.1/test/spectr.dat.0')
# V_g1=spc.Value('GAUGE',[1])
# V_g2=spc.Value('GAUGE',[2])
# V_mt1=spc.Value('MASS',[1000006])
# V_mt2=spc.Value('MASS',[2000006])
# V_mt=spc.Value('SMINPUTS',[6])
# V_At=spc.Value('EXTPAR',[11])
# V_MHU2=spc.Value('MSOFT',[22])
# V_MHD2=spc.Value('MSOFT',[21])
# V_MS2=spc.Value('NMSSMRUN',[10])
# V_TB=spc.Value('MINPAR',[3])
# V_Mz=spc.Value('MASS',[23])
# V_MU=spc.Value('NMSSMRUN',[5])
# V_AL=spc.Value('NMSSMRUN',[3])
# V_AK=spc.Value('NMSSMRUN',[4])
# V_L=spc.Value('NMSSMRUN',[1])
# V_K=spc.Value('NMSSMRUN',[2])
# V_mh2=spc.Value('MASS',[25])**2
# V_yt=2**0.5*V_mt/(sympy.sin(sympy.atan(V_TB))*246)
# #V_yt=V_mt/(sympy.sin(sympy.atan(V_TB))*174)
# # print('获取数值完成')


#Dpi=Delta_pi
DMHU2=sympy.Symbol('DMHU2')
DMHD2=sympy.Symbol('DMHD2')
DMS2=sympy.Symbol('DMS2')
DAL=sympy.Symbol('DAL')
DAK=sympy.Symbol('DAK')
DL=sympy.Symbol('DL')
DK=sympy.Symbol('DK')
Dyt=sympy.Symbol('Dyt')
#minimisation equation
E1=MHU2+MU**2+L**2*2*Mz2/(G*(1+TB2))+Mz2*(TB2-1)/(2*(1+TB2))-MU/(TB*L)*(AL*L+MU*K)+Wt
E2=MHD2+MU**2+L**2*2*Mz2*TB2/(G*(1+TB2))+Mz2*(1-TB2)/(2*(1+TB2))-MU*TB/L*(AL*L+MU*K)
E3=MS2+K*AK*MU/L+2*K**2*MU**2/L**2+2*L**2*Mz2/G-4*L*K*Mz2*TB/(G*(1+TB2))-2*L**2*AL*Mz2*TB/(MU*G*(1+TB2))
#matrix elements square
M11=Mz2*TB2/(1+TB2)+MU/(TB*L)*(AL*L+MU*K)+2*Wt
M22=Mz2/(1+TB2)+MU*TB*(AL*L+MU*K)/L
M33=2*Mz2*TB*L**2*AL/(MU*G*(1+TB2))+MU*K/L**2*(AK*L+4*MU*K)
M12=-MU/L*(AL*L+MU*K)+2*Mz2*TB/(G*(1+TB2))*(2*L**2-G/2)
M21=M12
M13=Mz*(2*L*MU*TB-(AL*L+2*MU*K))*sympy.root(2/(G*(1+TB2)),2)
M31=M13
M23=Mz*(2*L*MU-TB*(AL*L+2*MU*K))*sympy.root(2/(G*(1+TB2)),2)
M32=M23

#diagonalization of matrix  det_f
f=-(mh2)**3+(mh2)**2*(M11+M22+M33)-(mh2)*(M11*M22+M11*M33+M22*M33-M12*M21-M13*M31-M23*M32)+M11*M22*M33+M12*M23*M31+M13*M21*M32-M22*M13*M31-M12*M21*M33-M23*M32*M11
# print('质量平方矩阵对角化完成')

#Coefficient matrix
A11=sympy.diff(E1,Mz2)
A12=sympy.diff(E1,TB)
A13=sympy.diff(E1,MU)
A21=sympy.diff(E2,Mz2)
A22=sympy.diff(E2,TB)
A23=sympy.diff(E2,MU)
A31=sympy.diff(E3,Mz2)
A32=sympy.diff(E3,TB)
A33=sympy.diff(E3,MU)
#Augmented matrix term
D11=sympy.diff(E1,MHU2)*DMHU2+sympy.diff(E1,MHD2)*DMHD2+sympy.diff(E1,MS2)*DMS2+sympy.diff(E1,AL)*DAL+sympy.diff(E1,AK)*DAK+sympy.diff(E1,L)*DL+sympy.diff(E1,K)*DK+sympy.diff(E1,yt)*Dyt
D12=sympy.diff(E2,MHU2)*DMHU2+sympy.diff(E2,MHD2)*DMHD2+sympy.diff(E2,MS2)*DMS2+sympy.diff(E2,AL)*DAL+sympy.diff(E2,AK)*DAK+sympy.diff(E2,L)*DL+sympy.diff(E2,K)*DK+sympy.diff(E2,yt)*Dyt
D13=sympy.diff(E3,MHU2)*DMHU2+sympy.diff(E3,MHD2)*DMHD2+sympy.diff(E3,MS2)*DMS2+sympy.diff(E3,AL)*DAL+sympy.diff(E3,AK)*DAK+sympy.diff(E3,L)*DL+sympy.diff(E3,K)*DK+sympy.diff(E3,yt)*Dyt
#solve equations use Cramer's Rule
#det(A) det(D)
A=sympy.Matrix([[A11,A12,A13],[A21,A22,A23],[A31,A32,A33]])
D1=sympy.Matrix([[-D11,A12,A13],[-D12,A22,A23],[-D13,A32,A33]])
D2=sympy.Matrix([[A11,-D11,A13],[A21,-D12,A23],[A31,-D13,A33]])
D3=sympy.Matrix([[A11,A12,-D11],[A21,A22,-D12],[A31,A32,-D13]])
DET_A=A.det()
DET_D1=D1.det()
DET_D2=D2.det()
DET_D3=D3.det()
#DMz=Delta_Mz
DMz=DET_D1/DET_A
#DTB=Delta_TB
DTB=DET_D2/DET_A
#DMU=Delta_MU
DMU=DET_D3/DET_A
# print('Delta_Mz,TB,MU 计算完成')
#DMzi=d(DMz/Dpi)
DMz1=sympy.diff(DMz,DMHU2)
DMz2=sympy.diff(DMz,DMHD2)
DMz3=sympy.diff(DMz,DMS2)
DMz4=sympy.diff(DMz,DAL)
DMz5=sympy.diff(DMz,DAK)
DMz6=sympy.diff(DMz,DL)
DMz7=sympy.diff(DMz,DK)
DMz8=sympy.diff(DMz,Dyt)
# print('Delta_Mz_i 计算完成')
#d(DTB/Dpi)
DTB1=sympy.diff(DTB,DMHU2)
DTB2=sympy.diff(DTB,DMHD2)
DTB3=sympy.diff(DTB,DMS2)
DTB4=sympy.diff(DTB,DAL)
DTB5=sympy.diff(DTB,DAK)
DTB6=sympy.diff(DTB,DL)
DTB7=sympy.diff(DTB,DK)
DTB8=sympy.diff(DTB,Dyt)
# print('Delta_TanB_i 计算完成')
#d(DMU/Dpi)
DMU1=sympy.diff(DMU,DMHU2)
DMU2=sympy.diff(DMU,DMHD2)
DMU3=sympy.diff(DMU,DMS2)
DMU4=sympy.diff(DMU,DAL)
DMU5=sympy.diff(DMU,DAK)
DMU6=sympy.diff(DMU,DL)
DMU7=sympy.diff(DMU,DK)
DMU8=sympy.diff(DMU,Dyt)
# print('Delta_MU_i 计算完成')
#d(f/Mz)
dfMz=sympy.diff(f,Mz)
#d(f/TB)
dfTB=sympy.diff(f,TB)
#d(f/MU)
dfMU=sympy.diff(f,MU)
#d(f/mh2)
dfmh2=sympy.diff(f,mh2)
#d(f/pi)
dfMHU2=sympy.diff(f,MHU2)
dfMHD2=sympy.diff(f,MHD2)
dfMS2=sympy.diff(f,MS2)
dfAL=sympy.diff(f,AL)
dfAK=sympy.diff(f,AK)
dfL=sympy.diff(f,L)
dfK=sympy.diff(f,K)
dfyt=sympy.diff(f,yt)
dfpi=dfMHU2+dfMHD2+dfMS2+dfAL+dfAK+dfL+dfK+dfyt
# print('d(f/(Mz,TB,MU,mh2,pi)) 计算完成')

V_dic={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt,Mz2:V_Mz2}
# D_mh_FTi=pi/mh2*(dfMz*DMZi+dfTB*DTBi+dfMU*DMUi+dfpi)/dfmh2
# print('正在计算D_mh Fine Tuning...')
# mh_FT1=-MHU2/mh2*(dfMz*DMz1+dfTB*DTB1+dfMU*DMU1+dfpi)/dfmh2
# mh_FT2=-MHD2/mh2*(dfMz*DMz2+dfTB*DTB2+dfMU*DMU2+dfpi)/dfmh2
# mh_FT3=-MS2/mh2*(dfMz*DMz3+dfTB*DTB3+dfMU*DMU3+dfpi)/dfmh2
# mh_FT4=-AL/mh2*(dfMz*DMz4+dfTB*DTB4+dfMU*DMU4+dfpi)/dfmh2
# mh_FT5=-AK/mh2*(dfMz*DMz5+dfTB*DTB5+dfMU*DMU5+dfpi)/dfmh2
# mh_FT6=-L/mh2*(dfMz*DMz6+dfTB*DTB6+dfMU*DMU6+dfpi)/dfmh2
# mh_FT7=-K/mh2*(dfMz*DMz7+dfTB*DTB7+dfMU*DMU7+dfpi)/dfmh2
# mh_FT8=-yt/mh2*(dfMz*DMz8+dfTB*DTB8+dfMU*DMU8+dfpi)/dfmh2
# mh_V_FT1=mh_FT1.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT1:',mh_V_FT1)
# mh_V_FT2=mh_FT2.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT2:',mh_V_FT2)
# mh_V_FT3=mh_FT3.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT3:',mh_V_FT3)
# mh_V_FT4=mh_FT4.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT4:',mh_V_FT4)
# mh_V_FT5=mh_FT5.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT5:',mh_V_FT5)
# mh_V_FT6=mh_FT6.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT6:',mh_V_FT6)
# mh_V_FT7=mh_FT7.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT7:',mh_V_FT7)
# mh_V_FT8=mh_FT8.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
# print('mh_FT8:',mh_V_FT8)

# #D_Mz_FTi=pi/Mz*
Mz_FT1=MHU2/Mz2*DMz1
Mz_V_FT1=Mz_FT1.evalf(subs=V_dic)
print('Mz_FT1:',Mz_V_FT1)

Mz_FT2=MHD2/Mz2*DMz2
Mz_V_FT2=Mz_FT2.evalf(subs=V_dic)
print('Mz_FT2:',Mz_V_FT2)

Mz_FT3=MS2/Mz2*DMz3
Mz_V_FT3=Mz_FT3.evalf(subs=V_dic)
print('Mz_FT3:',Mz_V_FT3)

Mz_FT4=AL/Mz2*DMz4
Mz_V_FT4=Mz_FT4.evalf(subs=V_dic)
print('Mz_FT4:',Mz_V_FT4)

Mz_FT5=AK/Mz2*DMz5
Mz_V_FT5=Mz_FT5.evalf(subs=V_dic)
print('Mz_FT5:',Mz_V_FT5)

Mz_FT6=L/Mz2*DMz6
Mz_V_FT6=Mz_FT6.evalf(subs=V_dic)
print('Mz_FT6:',Mz_V_FT6)

Mz_FT7=K/Mz2*DMz7
Mz_V_FT7=Mz_FT7.evalf(subs=V_dic)
print('Mz_FT7:',Mz_V_FT7)

Mz_FT8=yt/Mz2*DMz8
Mz_V_FT8=Mz_FT8.evalf(subs=V_dic)
print('Mz_FT8:',Mz_V_FT8)

# print('DTB1:',DTB1.evalf(subs=V_dic))
# print('DTB2:',DTB2.evalf(subs=V_dic))
# print('DTB3:',DTB3.evalf(subs=V_dic))
# print('DTB4:',DTB5.evalf(subs=V_dic))
# print('DTB5:',DTB8.evalf(subs=V_dic))
# print('DTB6:',DTB4.evalf(subs=V_dic))
# print('DTB7:',DTB6.evalf(subs=V_dic))
# print('DTB8:',DTB7.evalf(subs=V_dic))
# print('DMU1:',DMU1.evalf(subs=V_dic))
# print('DMU2:',DMU2.evalf(subs=V_dic))
# print('DMU3:',DMU3.evalf(subs=V_dic))
# print('DMU4:',DMU5.evalf(subs=V_dic))
# print('DMU5:',DMU8.evalf(subs=V_dic))
# print('DMU6:',DMU4.evalf(subs=V_dic))
# print('DMU7:',DMU6.evalf(subs=V_dic))
# print('DMU8:',DMU7.evalf(subs=V_dic))
# print('M11:',M11.evalf(subs=V_dic))
# print('M22:',M22.evalf(subs=V_dic))
# print('M33:',M33.evalf(subs=V_dic))
# print('M12:',M12.evalf(subs=V_dic))
# print('M13:',M13.evalf(subs=V_dic))
# print('M23:',M23.evalf(subs=V_dic))
# print('dfmh2:',dfmh2.evalf(subs=V_dic))
# print('dfMz:',dfMz.evalf(subs=V_dic))
# print('dfTB:',dfTB.evalf(subs=V_dic))
# print('dfMU:',dfMU.evalf(subs=V_dic))
# print('dfMHU2:',dfMHU2.evalf(subs=V_dic))
# print('dfMHD2:',dfMHD2.evalf(subs=V_dic))
# print('dfMS2:',dfMS2.evalf(subs=V_dic))
# print('dfAL:',dfAL.evalf(subs=V_dic))
# print('dfAK:',dfAK.evalf(subs=V_dic))
# print('dfL:',dfL.evalf(subs=V_dic))
# print('dfK:',dfK.evalf(subs=V_dic))
# print('dfyt:',dfyt.evalf(subs=V_dic))
# print('dfpi:',dfpi.evalf(subs=V_dic))

# print('dM11Z:',sympy.diff(M11,Mz).evalf(subs=V_dic))
# print('dM22Z:',sympy.diff(M22,Mz).evalf(subs=V_dic))
# print('dM33Z:',sympy.diff(M33,Mz).evalf(subs=V_dic))
# print('dM12Z:',sympy.diff(M12,Mz).evalf(subs=V_dic))
# print('dM13Z:',sympy.diff(M13,Mz).evalf(subs=V_dic))
# print('dM23Z:',sympy.diff(M23,Mz).evalf(subs=V_dic))
# print('dM11TB:',sympy.diff(M11,TB).evalf(subs=V_dic))
# print('dM22TB:',sympy.diff(M22,TB).evalf(subs=V_dic))
# print('dM33TB:',sympy.diff(M33,TB).evalf(subs=V_dic))
# print('dM12TB:',sympy.diff(M12,TB).evalf(subs=V_dic))
# print('dM13TB:',sympy.diff(M13,TB).evalf(subs=V_dic))
# print('dM23TB:',sympy.diff(M23,TB).evalf(subs=V_dic))
# print('dM11MU:',sympy.diff(M11,MU).evalf(subs=V_dic))
# print('dM22MU:',sympy.diff(M22,MU).evalf(subs=V_dic))
# print('dM33MU:',sympy.diff(M33,MU).evalf(subs=V_dic))
# print('dM12MU:',sympy.diff(M12,MU).evalf(subs=V_dic))
# print('dM13MU:',sympy.diff(M13,MU).evalf(subs=V_dic))
# print('dM23MU:',sympy.diff(M23,MU).evalf(subs=V_dic))
# print('dM11L;',sympy.diff(M11,L).evalf(subs=V_dic))
# print('dM22L;',sympy.diff(M22,L).evalf(subs=V_dic))
# print('dM33L;',sympy.diff(M33,L).evalf(subs=V_dic))
# print('dM12L;',sympy.diff(M12,L).evalf(subs=V_dic))
# print('dM13L;',sympy.diff(M13,L).evalf(subs=V_dic))
# print('dM23L;',sympy.diff(M23,L).evalf(subs=V_dic))
# print('dM11K;',sympy.diff(M11,K).evalf(subs=V_dic))
# print('dM22K;',sympy.diff(M22,K).evalf(subs=V_dic))
# print('dM33K;',sympy.diff(M33,K).evalf(subs=V_dic))
# print('dM12K;',sympy.diff(M12,K).evalf(subs=V_dic))
# print('dM13K;',sympy.diff(M13,K).evalf(subs=V_dic))
# print('dM23K;',sympy.diff(M23,K).evalf(subs=V_dic))
# print('dM11AL:',sympy.diff(M11,AL).evalf(subs=V_dic))
# print('dM22AL:',sympy.diff(M22,AL).evalf(subs=V_dic))
# print('dM33AL:',sympy.diff(M33,AL).evalf(subs=V_dic))
# print('dM12AL:',sympy.diff(M12,AL).evalf(subs=V_dic))
# print('dM13AL:',sympy.diff(M13,AL).evalf(subs=V_dic))
# print('dM23AL:',sympy.diff(M23,AL).evalf(subs=V_dic))
# print('dM11AK:',sympy.diff(M11,AK).evalf(subs=V_dic))
# print('dM22AK:',sympy.diff(M22,AK).evalf(subs=V_dic))
# print('dM33AK:',sympy.diff(M33,AK).evalf(subs=V_dic))
# print('dM12AK:',sympy.diff(M12,AK).evalf(subs=V_dic))
# print('dM13AK:',sympy.diff(M13,AK).evalf(subs=V_dic))
# print('dM23AK:',sympy.diff(M23,AK).evalf(subs=V_dic))
# print('dM11yt:',sympy.diff(M11,yt).evalf(subs=V_dic))
# print('dM22yt:',sympy.diff(M22,yt).evalf(subs=V_dic))
# print('dM33yt:',sympy.diff(M33,yt).evalf(subs=V_dic))
# print('dM12yt:',sympy.diff(M12,yt).evalf(subs=V_dic))
# print('dM13yt:',sympy.diff(M13,yt).evalf(subs=V_dic))
# print('dM23yt:',sympy.diff(M23,yt).evalf(subs=V_dic))


# f1=(mh2)**2*(M11+M22+M33)
# f2=-(mh2)*(M11*M22+M11*M33+M22*M33-M12*M21-M13*M31-M23*M32)
# f3=M11*M22*M33+M12*M23*M31+M13*M21*M32-M22*M13*M31-M12*M21*M33-M23*M32*M11

# df1Mz=sympy.diff(f1,Mz)
# df2Mz=sympy.diff(f2,Mz)
# df3Mz=sympy.diff(f3,Mz)

# print('df1Mz:',df1Mz.evalf(subs=V_dic))
# print('df2Mz:',df2Mz.evalf(subs=V_dic))
# print('df3Mz:',df3Mz.evalf(subs=V_dic))

# f31=M11*M22*M33
# f32=M12*M23*M31+M13*M21*M32
# f33=-M22*M13*M31-M12*M21*M33-M23*M32*M11

# df31Mz=sympy.diff(f31,Mz)
# df32Mz=sympy.diff(f32,Mz)
# df33Mz=sympy.diff(f33,Mz)

# print('df31Mz:',df31Mz.evalf(subs=V_dic))
# print('df32Mz:',df32Mz.evalf(subs=V_dic))
# print('df33Mz:',df33Mz.evalf(subs=V_dic))

# df3_Mz=df31Mz+df32Mz+df33Mz
# print('df3_Mz:',df3_Mz.evalf(subs=V_dic))

print('计算完毕')
end=time.clock()
print('用时',end-start)