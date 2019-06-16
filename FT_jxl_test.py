import time
start=time.clock()
#import xslha
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
MU=sympy.Symbol('MU')
AL=sympy.Symbol('AL')
AK=sympy.Symbol('AK')
L=sympy.Symbol('L')	
K=sympy.Symbol('K')	
yt=sympy.Symbol('yt')
origin=3*yt**4/(8*(sympy.pi)**2)*2*Mz**2*TB2/(G*(1+TB2))
#Wt=origin*sympy.log(mt1**2/mt2**2)
Wt=origin*(A+B*(At+MU/TB)+C*(At+MU/TB)**2)
print('导入参数完成')
#获得数值
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
print('获取数值完成')
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
E1=MHU2+MU**2+L**2*2*Mz**2/(G*(1+TB2))+Mz**2*(TB2-1)/(2*(1+TB2))-MU/(TB*L)*(AL*L+MU*K)+Wt
E2=MHD2+MU**2+L**2*2*Mz**2*TB2/(G*(1+TB2))+Mz**2*(1-TB2)/(2*(1+TB2))-MU*TB/L*(AL*L+MU*K)
E3=MS2+K*AK*MU/L+2*K**2*MU**2/L**2+2*L**2*Mz**2/G-4*L*K*Mz**2*TB/(G*(1+TB2))-2*L**2*AL*Mz**2*TB/(MU*G*(1+TB2))
#matrix elements square
M11=Mz**2*TB2/(1+TB2)+MU/(TB*L)*(AL*L+MU*K)+Wt
M22=Mz**2/(1+TB2)+MU*TB*(AL*L+MU*K)/L
M33=2*Mz**2*TB*L**2*AL/(MU*G*(1+TB2))+MU*K/L**2*(AK*L+4*MU*K)
M12=-MU/L*(AL*L+MU*K)+2*Mz**2*TB/(G*(1+TB2))*(2*L**2-G/2)
M21=M12
M13=Mz*(2*L*MU*TB-(AL*L+2*MU*K))*sympy.root(2/(G*(1+TB2)),2)
M31=M13
M23=Mz*(2*L*MU-TB*(AL*L+2*MU*K))*sympy.root(2/(G*(1+TB2)),2)
M32=M23

#diagonalization of matrix  det_f
f=-(mh2)**3+(mh2)**2*(M11+M22+M33)-(mh2)*(M11*M22+M11*M33+M22*M33-M12*M21-M13*M31-M23*M32)+M11*M22*M33+M12*M23*M31+M13*M21*M32-M22*M13*M31-M12*M21*M33-M23*M32*M11
print('质量平方矩阵对角化完成')
V_mh2=15703.778547419524
#Coefficient matrix
A11=sympy.diff(E1,Mz)
A12=sympy.diff(E1,TB)
A13=sympy.diff(E1,MU)
A21=sympy.diff(E2,Mz)
A22=sympy.diff(E2,TB)
A23=sympy.diff(E2,MU)
A31=sympy.diff(E3,Mz)
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
print('Delta_Mz,TB,MU 计算完成')
#DMzi=d(DMz/Dpi)
DMz1=sympy.diff(DMz,DMHU2)
DMz2=sympy.diff(DMz,DMHD2)
DMz3=sympy.diff(DMz,DMS2)
DMz4=sympy.diff(DMz,DAL)
DMz5=sympy.diff(DMz,DAK)
DMz6=sympy.diff(DMz,DL)
DMz7=sympy.diff(DMz,DK)
DMz8=sympy.diff(DMz,Dyt)
print('Delta_Mz_i 计算完成')
#d(DTB/Dpi)
DTB1=sympy.diff(DTB,DMHU2)
DTB2=sympy.diff(DTB,DMHD2)
DTB3=sympy.diff(DTB,DMS2)
DTB4=sympy.diff(DTB,DAL)
DTB5=sympy.diff(DTB,DAK)
DTB6=sympy.diff(DTB,DL)
DTB7=sympy.diff(DTB,DK)
DTB8=sympy.diff(DTB,Dyt)
print('Delta_TanB_i 计算完成')
#d(DMU/Dpi)
DMU1=sympy.diff(DMU,DMHU2)
DMU2=sympy.diff(DMU,DMHD2)
DMU3=sympy.diff(DMU,DMS2)
DMU4=sympy.diff(DMU,DAL)
DMU5=sympy.diff(DMU,DAK)
DMU6=sympy.diff(DMU,DL)
DMU7=sympy.diff(DMU,DK)
DMU8=sympy.diff(DMU,Dyt)
print('Delta_MU_i 计算完成')
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
print('d(f/(Mz,TB,MU,mh2,pi)) 计算完成')

#D_mh_FTi=pi/mh2*(dfMz*DMZi+dfTB*DTBi+dfMU*DMUi+dfpi)/dfmh2
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
Mz_FT1=MHU2/Mz*DMz1
Mz_V_FT1=Mz_FT1.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT1:',Mz_V_FT1)

Mz_FT2=MHD2/Mz*DMz2
Mz_V_FT2=Mz_FT2.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT2:',Mz_V_FT2)

Mz_FT3=MS2/Mz*DMz3
Mz_V_FT3=Mz_FT3.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT3:',Mz_V_FT3)

Mz_FT4=AL/Mz*DMz4
Mz_V_FT4=Mz_FT4.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT4:',Mz_V_FT4)

Mz_FT5=AK/Mz*DMz5
Mz_V_FT5=Mz_FT5.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT5:',Mz_V_FT5)

Mz_FT6=L/Mz*DMz6
Mz_V_FT6=Mz_FT6.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT6:',Mz_V_FT6)

Mz_FT7=K/Mz*DMz7
Mz_V_FT7=Mz_FT7.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT7:',Mz_V_FT7)

Mz_FT8=yt/Mz*DMz8
Mz_V_FT8=Mz_FT8.evalf(subs={mh2:V_mh2,g1:V_g1,g2:V_g2,mt1:V_mt1,mt2:V_mt2,mt:V_mt,At:V_At,MHU2:V_MHU2,MHD2:V_MHD2,MS2:V_MS2,TB:V_TB,Mz:V_Mz,MU:V_MU,AL:V_AL,AK:V_AK,L:V_L,K:V_K,yt:V_yt})
print('Mz_FT8:',Mz_V_FT8)

print('计算完毕')
end=time.clock()
print('用时',end-start)