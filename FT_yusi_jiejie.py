
import numpy as np
import math
import re

mh2=float(1.25314718E+02)**2
input_parameter={'At':2000.0,'mu':300.0,'tanbeta':float(1.14531792E+01),'Akappa':200.0,'Lambda':float(2.28888827E-01),'kappa':float(-1.39604166E-01),'Alambda':2000.0,'Yt':0.003,
'mz':float(9.11887000E+01),'mt':float(1.73500000E+02),'g1':float(3.61496259E-01),'g2':float(6.71344778E-01),'mhd2':float(3.45021908E+06),'mhu2':float(-3.40835683E+05),'ms2':float(1.79817995E+04),'mst1':float(2.01037642E+03),'mst2':float(2.12820584E+03),'mh2':15703.778547419524}


def ZH_Finetuning(input_parameter):

    # global mz,mt,g1,g2,mhd2,mhdu2,ms2,mst1,mst2,mh2,At,mu,tb,ak,l,k,al,Yt,g2,lh,kh,b,tb2,ct,a1,a2,a3,c

    mz  = input_parameter['mz']
    mt  = input_parameter['mt']
    g1  = input_parameter['g1']
    g2  = input_parameter['g2']
    mhd2  = input_parameter['mhd2']
    mhu2  = input_parameter['mhu2']
    ms2  = input_parameter['ms2']
    mst1  = input_parameter['mst1']
    mst2  = input_parameter['mst2']
    mh2  = input_parameter['mh2']
    At  = input_parameter['At']
    mu  = input_parameter['mu']
    tb  = input_parameter['tanbeta']
    ak  = input_parameter['Akappa']
    l   = input_parameter['Lambda']
    k   = input_parameter['kappa']
    al  = input_parameter['Alambda']
    Yt  = input_parameter['Yt']



    # global g2,tb,l,k,al,ak,lh,kh,b,tb2,ct,a1,a2,a3,c
    g2  = (g1**2+g2**2)/2
    lh  = l**2/g2
    kh  = k/l
    b   = al+kh*mu
    tb2 = tb**2

    c   = (3*Yt**4)/(8*math.pi**2*g2)
    
    a1  = np.log(mst1*mst2/mt**2)
    a2  = At/(mst1**2-mst2**2)*np.log(mst1**2/mst2**2)
    a3  = At**2/(mst1**2-mst2**2)**2*(1-(mst1**2+mst2**2)/(mst1**2-mst2**2)*np.log(mst1/mst2))
    ct  = c*(a1+a2*(At+mu/tb)+a3*(At+mu/tb)**2)

    def mz_derivative():
        # DiT indicates Ei derivativing of tanbeta. DiT = d(i)/d(tb)
        D1T= 2*mz**2*tb/(1+tb2)**2*(1-lh+ct)+(mu*b)/tb2+c*mz**2*tb2/(1+tb2)*(-a2*mu/tb2-2*a3*mu/tb2*(At+mu/tb))
        D2T= 2*mz**2*tb/(1+tb2)**2*(lh-1)-mu*b
        D3T= lh*mz**2*(tb2-1)/(1+tb2)**2*(b+kh*mu)/mu

        # DiM indicates Ei derivativing of mu. DiM = d(i)/d(mu)
        D1M= 2*mu-(b+kh*mu)/tb+c*mz**2*tb2/(1+tb2)*(a2/tb+2*a3/tb*(At+mu/tb))
        D2M= 2*mu-(b+kh*mu)*tb
        D3M= kh*ak+4*kh**2*mu+lh*mz**2*tb/(1+tb2)*al/mu**2

        # DiZ indicates Ei derivativing of mz. DiZ = d(i)/d(mz)
        D1Z= 2*mz*(lh+(tb2-1)/2+ct*tb2)/(1+tb2)
        D2Z= 2*mz*(lh*tb2+(1-tb2)/2)/(1+tb2)
        D3Z= 2*mz*lh*(1-tb/(1+tb2)*(b+kh*mu)/mu)

        #Ai=d(j)/d(tb)*d(k)/d(mu)-d(k)/d(tb)*d(j)/d(mu)
        A1=-D2T*D3M+D2M*D3T
        A2=-D3T*D1M+D3M*D1T
        A3=-D1T*D2M+D1M*D2T

        #Bi=d(j)/d(mu)*d(k)/d(mz)-d(k)/d(mu)*d(j)/d(mz)
        B1=-D2M*D3Z+D2Z*D3M
        B2=-D3M*D1Z+D3Z*D1M
        B3=-D1M*D2Z+D1Z*D2M

        #Ci=d(j)/d(tb)*d(k)/d(mz)-d(k)/d(tb)*d(j)/d(mz)
        C1=-(-D2T*D3Z+D2Z*D3T)
        C2=-(-D3T*D1Z+D3Z*D1T)
        C3=-(-D1T*D2Z+D1Z*D2T)

        #Determinant of the linear system
        D=-(A1*D1Z+A2*D2Z+A3*D3Z)

        #dmz=d(mz)/d(pi)
        #pi={mhu2,mhd2,ms2,Akappa,Yt,Alambda,Lambda,kappa}
        dmz1=A1/D                            ###dmz1 indicates d(mz)/d(mhu2)
        dmz2=A2/D                            ###dmz2 indicates d(mz)/d(mhd2)
        dmz3=A3/D                            ###dmz3 indicates d(mz)/d(ms2)
        dmz4=A3/D*mu*kh                      ###dmz4 indicates d(mz)/d(Akappa)
        dmz5=A1/(D*Yt)*4*ct*mz**2*tb2/(1+tb2)        ###dmz5 indicates d(mz)/d(Yt)
        dmz6=-(A1*mu/tb+A2*mu*tb+A3*(mz**2*tb*lh)/(mu*(1+tb2)))/D       ###dmz6 indicates d(mz)/d(Alambda)
        dmz7=(A1*(2*mz**2*lh/(1+tb2)+kh*mu**2/tb)+A2*(2*mz**2*tb2*lh/(1+tb2)+kh*mu**2*tb)+A3*(2*mz**2*lh-2*lh*mz**2*tb/(1+tb2)*(kh+al/mu)-mu*kh*(ak+4*mu*kh)))/(D*l)     ###dmz7 indicates d(mz)/d(lambda)
        dmz8=-kh*(A1*(mu**2/tb)+A2*(mu**2*tb)+A3*(2*mz**2*tb*lh/(1+tb2)-mu*ak-4*mu**2*kh))/(D*k)    ###dmz8 indicates d(mz)/d(kappa)

        #dtb=d(tb)/d(pi)
        #pi={mhu2,mhd2,ms2,Akappa,Yt,Alambda,Lambda,kappa}
        dtb1=B1/D                            ###dtb1 indicates d(tb)/d(mhu2)
        dtb2=B2/D                            ###dtb2 indicates d(tb)/d(mhd2)
        dtb3=B3/D                            ###dtb3 indicates d(tb)/d(ms2)
        dtb4=B3/D*mu*kh                      ###dtb4 indicates d(tb)/d(Akappa)
        dtb5=B1/(D*Yt)*4*ct*mz**2*tb2/(1+tb2)        ###dtb5 indicates d(tb)/d(Yt)
        dtb6=-(B1*mu/tb+B2*mu*tb+B3*(mz**2*tb*lh)/(mu*(1+tb2)))/D       ###dtb6 indicates d(tb)/d(Alambda)
        dtb7=(B1*(2*mz**2*lh/(1+tb2)+kh*mu**2/tb)+B2*(2*mz**2*tb2*lh/(1+tb2)+kh*mu**2*tb)+B3*(2*mz**2*lh-2*lh*mz**2*tb/(1+tb2)*(kh+al/mu)-mu*kh*(ak+4*mu*kh)))/(D*l)     ###dtb7 indicates d(tb)/d(lambda)
        dtb8=-kh*(B1*(mu**2/tb)+B2*(mu**2*tb)+B3*(2*mz**2*tb*lh/(1+tb2)-mu*ak-4*mu**2*kh))/(D*k)    ###dtb8 indicates d(tb)/d(kappa)

        #dmu=d(mu)/d(pi)
        #pi={mhu2,mhd2,ms2,Akappa,Yt,Alambda,Lambda,kappa}
        dmu1=C1/D                            ###dmu1 indicates d(mu)/d(mhu2)
        dmu2=C2/D                            ###dmu2 indicates d(mu)/d(mhd2)
        dmu3=C3/D                            ###dmu3 indicates d(mu)/d(ms2)
        dmu4=C3/D*mu*kh                      ###dmu4 indicates d(mu)/d(Akappa)
        dmu5=C1/(D*Yt)*4*ct*mz**2*tb2/(1+tb2)        ###dmu5 indicates d(mu)/d(Yt)
        dmu6=-(C1*mu/tb+C2*mu*tb+C3*(mz**2*tb*lh)/(mu*(1+tb2)))/D       ###dmu6 indicates d(mu)/d(Alambda)
        dmu7=(C1*(2*mz**2*lh/(1+tb2)+kh*mu**2/tb)+C2*(2*mz**2*tb2*lh/(1+tb2)+kh*mu**2*tb)+C3*(2*mz**2*lh-2*lh*mz**2*tb/(1+tb2)*(kh+al/mu)-mu*kh*(ak+4*mu*kh)))/(D*l)     ###dmu7 indicates d(mu)/d(lambda)
        dmu8=-kh*(C1*(mu**2/tb)+C2*(mu**2*tb)+C3*(2*mz**2*tb*lh/(1+tb2)-mu*ak-4*mu**2*kh))/(D*k)    ###dmu8 indicates d(mu)/d(kappa)
        
        FT_mz1=np.abs(mhu2/mz*dmz1)
        FT_mz2=np.abs(mhd2/mz*dmz2)
        FT_mz3=np.abs(ms2/mz*dmz3)
        FT_mz4=np.abs(ak/mz*dmz4)
        FT_mz5=np.abs(Yt/mz*dmz5)
        FT_mz6=np.abs(al/mz*dmz6)
        FT_mz7=np.abs(l/mz*dmz7)
        FT_mz8=np.abs(k/mz*dmz8)

        FT_mz={'FT_mz1':FT_mz1,'FT_mz2':FT_mz2,'FT_mz3':FT_mz3,'FT_mz4':FT_mz4,
        'FT_mz5':FT_mz5,'FT_mz6':FT_mz6,'FT_mz7':FT_mz7,'FT_mz8':FT_mz8}

        delta_mz={'dmz1':dmz1,'dmz2':dmz2,'dmz3':dmz3,'dmz4':dmz4,'dmz5':dmz5,'dmz6':dmz6,'dmz7':dmz7,'dmz8':dmz8}

        delta_mz_all={'dmz1':dmz1,'dmz2':dmz2,'dmz3':dmz3,'dmz4':dmz4,'dmz5':dmz5,'dmz6':dmz6,'dmz7':dmz7,'dmz8':dmz8,
        'dtb1':dtb1,'dtb2':dtb2,'dtb3':dtb3,'dtb4':dtb4,'dtb5':dtb5,'dtb6':dtb6,'dtb7':dtb7,'dtb8':dtb8,
        'dmu1':dmu1,'dmu2':dmu2,'dmu3':dmu3,'dmu4':dmu4,'dmu5':dmu5,'dmu6':dmu6,'dmu7':dmu7,'dmu8':dmu8}
    
        # return FT_mz,delta_mz,delta_mz_all

        def mh_derivative():
            '''
            E=-(mh2)**3+(mh2)**2*(M11+M22+M33)-mh2*(M11*M22+M11*M33+M22*M33-M12*M21-M13*M31-M23*M32)
            +M11*M22*M33+2*M12*M13*M23-M22*M13*M31-M12*M21*M33-M23*M32*M11
            '''

            ### mass matrix   basis is (hu,hd,s)
            M11 = mz**2*tb2/(1+tb2)*(1+2*ct)+mu*b/tb
            M22 = mz**2/(1+tb2)+mu*tb*b
            M33 = lh*mz**2*al*tb/(mu*(1+tb2))+kh*mu*(ak+4*kh*mu)
            M12 = -mu*b+mz**2*tb*(2*lh-1)/(1+tb2)
            M13 = mz*(2*mu*tb-b-kh*mu)*math.sqrt(lh/(1+tb2))
            M23 = mz*(2*mu-tb*(b+kh*mu))*math.sqrt(lh/(1+tb2))

            ###dM=d(M)/d(pi)
            #pi={mhu2,mhd2,ms2,Akappa,Yt,Alambda,Lambda,kappa}
            #i=1 indicates mhu2   ,mhd2,ms2
            dM11_1  = 0
            dM22_1  = 0
            dM33_1  = 0
            dM12_1  = 0
            dM13_1  = 0
            dM23_1  = 0

            #i=2 indicates mhd2
            dM11_2  = 0
            dM22_2  = 0
            dM33_2  = 0
            dM12_2  = 0
            dM13_2  = 0
            dM23_2  = 0

            #i=3 indicates ms2
            dM11_3  = 0
            dM22_3  = 0
            dM33_3  = 0
            dM12_3  = 0
            dM13_3  = 0
            dM23_3  = 0
            
            #i=4 indicates Akappa
            dM11_4  = 0
            dM22_4  = 0
            dM33_4  = mu*kh
            dM12_4  = 0
            dM13_4  = 0
            dM23_4  = 0

            #i=5 indicates Yt
            dM11_5  = 8*ct*mz**2*tb2/((1+tb2)*Yt)
            dM22_5  = 0
            dM33_5  = 0
            dM12_5  = 0
            dM13_5  = 0
            dM23_5  = 0

            #i=6 indicates Alambda
            dM11_6  = mu/tb
            dM22_6  = mu*tb
            dM33_6  = mz**2*tb*lh/(mu*(1+tb2))
            dM12_6  = -mu
            dM13_6  = -mz*math.sqrt(lh/(1+tb2))
            dM23_6  = -mz*tb*math.sqrt(lh/(1+tb2))

            #i=7 indicates Lambda
            dM11_7  = -mu**2*kh/(tb*l)
            dM22_7  = -mu**2*kh*tb/l
            dM33_7  = 2*mz**2*tb*al*lh/(mu*l*(1+tb2))-mu*kh/l*(ak+8*mu*kh)
            dM12_7  = mu**2*kh/l+4*mz**2*tb*lh/(l*(1+tb2))
            dM13_7  = mz*(2*mu*tb-al)/l*math.sqrt(lh/(1+tb2))
            dM23_7  = mz*(2*mu-tb*al)/l*math.sqrt(lh/(1+tb2))

            #i=8 indicates kappa
            dM11_8  = mu**2/(tb*l)
            dM22_8  = mu**2*tb/l
            dM33_8  = mu/l*(ak+8*mu*kh)
            dM12_8  = -mu**2/l
            dM13_8  = -2*mz*mu/l*math.sqrt(lh/(1+tb2))
            dM23_8  = -2*mz*mu*tb/l*math.sqrt(lh/(1+tb2))

            ###dM**z=d(M)/d(mz)
            dM11z=2*mz*tb2/(1+tb2)+4*mz*tb2/(1+tb2)*ct
            dM22z=2*mz/(1+tb2)
            dM33z=2*mz*tb*al*lh/(mu*(1+tb2))
            dM12z=2*mz*tb/(1+tb2)*(2*lh-1)
            dM13z=(2*mu*tb-(b+mu*kh))*math.sqrt(lh/(1+tb2))
            dM23z=(2*mu-tb*(b+mu*kh))*math.sqrt(lh/(1+tb2))

            ###dM**t=d(M)/d(tb)
            dM11t=2*mz**2*tb/(1+tb2)**2-mu/tb2*b+4*mz**2*tb*ct/(1+tb2)**2+2*c*mz**2*tb2/(1+tb2)*(-a2*mu/tb2-2*a3*mu/tb2*(At+mu/tb))
            dM22t=-2*tb*mz**2/(1+tb2)**2+mu*b
            dM33t=mz**2*(1-tb2)*al*lh/(mu*(1+tb2)**2)
            dM12t=mz**2*(1-tb2)*(2*lh-1)/(1+tb2)**2
            dM13t=mz*(2*mu+(b+mu*kh)*tb)/(1+tb2)*math.sqrt(lh/(1+tb2))
            dM23t=-mz*(2*mu*tb+(b+mu*kh))/(1+tb2)*math.sqrt(lh/(1+tb2))

            ###dM**m=d(M)/d(mu)
            dM11m=(b+mu*kh)/tb+2*c*mz**2*tb2/(1+tb2)*(a2/tb+2*a3*(At+mu/tb)/tb)
            dM22m=tb*(b+mu*kh)
            dM33m=-mz**2*tb*al*lh/(mu**2*(1+tb2))+kh*(ak+8*mu*kh)
            dM12m=-(b+mu*kh)
            dM13m=mz*(2*tb-2*kh)*math.sqrt(lh/(1+tb2))
            dM23m=mz*(2-2*tb*kh)*math.sqrt(lh/(1+tb2))

            ###dEh=d(E)/d(mh2)
            dEh=-3*mh2**2+2*mh2*(M11+M22+M33)-(M11*M22+M11*M33+M22*M33-M12**2-M13**2-M23**2)

            ###dEz=d(E)/d(mz)
            dEz=mh2**2*(dM11z+dM22z+dM33z)-mh2*((dM22z+dM33z)*M11+(dM11z+dM33z)*M22+(dM11z+dM22z)*M33
            -2*(dM12z*M12+dM13z*M13+dM23z*M23))+(dM11z*M22*M33+dM22z*M11*M33+dM33z*M11*M22)
            +2*(dM12z*M13*M23+dM13z*M12*M23+dM23z*M12*M13)-(dM11z*M23**2+dM22z*M13**2+dM33z*M12**2)
            -2*(dM12z*M12*M33+dM13z*M13*M22+dM23z*M23*M11)

            ###dEt=d(E)/d(tb)
            dEt=mh2**2*(dM11t+dM22t+dM33t)-mh2*((dM22t+dM33t)*M11+(dM11t+dM33t)*M22+(dM11t+dM22t)*M33
            -2*(dM12t*M12+dM13t*M13+dM23t*M23))+(dM11t*M22*M33+dM22t*M11*M33+dM33t*M11*M22)
            +2*(dM12t*M13*M23+dM13t*M12*M23+dM23t*M12*M13)-(dM11t*M23**2+dM22t*M13**2+dM33t*M12**2)
            -2*(dM12t*M12*M33+dM13t*M13*M22+dM23t*M23*M11)

            ###dEm=d(E)/d(mu)
            dEm=mh2**2*(dM11m+dM22m+dM33m)-mh2*((dM22m+dM33m)*M11+(dM11m+dM33m)*M22+(dM11m+dM22m)*M33
            -2*(dM12m*M12+dM13m*M13+dM23m*M23))+(dM11m*M22*M33+dM22m*M11*M33+dM33m*M11*M22)
            +2*(dM12m*M13*M23+dM13m*M12*M23+dM23m*M12*M13)-(dM11m*M23**2+dM22m*M13**2+dM33m*M12**2)
            -2*(dM12m*M12*M33+dM13m*M13*M22+dM23m*M23*M11)

            ###dEp=d(E)/d(pi)
            dEp_1=mh2**2*(dM11_1+dM22_1+dM33_1)-mh2*((dM22_1+dM33_1)*M11+(dM11_1+dM33_1)*M22+(dM11_1+dM22_1)*M33
            -2*(dM12_1*M12+dM13_1*M13+dM23_1*M23))+(dM11_1*M22*M33+dM22_1*M11*M33+dM33_1*M11*M22)
            +2*(dM12_1*M13*M23+dM13_1*M12*M23+dM23_1*M12*M13)-(dM11_1*M23**2+dM22_1*M13**2+dM33_1*M12**2)
            -2*(dM12_1*M12*M33+dM13_1*M13*M22+dM23_1*M23*M11)

            dEp_2=mh2**2*(dM11_2+dM22_2+dM33_2)-mh2*((dM22_2+dM33_2)*M11+(dM11_2+dM33_2)*M22+(dM11_2+dM22_2)*M33
            -2*(dM12_2*M12+dM13_2*M13+dM23_2*M23))+(dM11_2*M22*M33+dM22_2*M11*M33+dM33_2*M11*M22)
            +2*(dM12_2*M13*M23+dM13_2*M12*M23+dM23_2*M12*M13)-(dM11_2*M23**2+dM22_2*M13**2+dM33_2*M12**2)
            -2*(dM12_2*M12*M33+dM13_2*M13*M22+dM23_2*M23*M11)

            dEp_3=mh2**2*(dM11_3+dM22_3+dM33_3)-mh2*((dM22_3+dM33_3)*M11+(dM11_3+dM33_3)*M22+(dM11_3+dM22_3)*M33
            -2*(dM12_3*M12+dM13_3*M13+dM23_3*M23))+(dM11_3*M22*M33+dM22_3*M11*M33+dM33_3*M11*M22)
            +2*(dM12_3*M13*M23+dM13_3*M12*M23+dM23_3*M12*M13)-(dM11_3*M23**2+dM22_3*M13**2+dM33_3*M12**2)
            -2*(dM12_3*M12*M33+dM13_3*M13*M22+dM23_3*M23*M11)

            dEp_4=mh2**2*(dM11_4+dM22_4+dM33_4)-mh2*((dM22_4+dM33_4)*M11+(dM11_4+dM33_4)*M22+(dM11_4+dM22_4)*M33
            -2*(dM12_4*M12+dM13_4*M13+dM23_4*M23))+(dM11_4*M22*M33+dM22_4*M11*M33+dM33_4*M11*M22)
            +2*(dM12_4*M13*M23+dM13_4*M12*M23+dM23_4*M12*M13)-(dM11_4*M23**2+dM22_4*M13**2+dM33_4*M12**2)
            -2*(dM12_4*M12*M33+dM13_4*M13*M22+dM23_4*M23*M11)

            dEp_5=mh2**2*(dM11_5+dM22_5+dM33_5)-mh2*((dM22_5+dM33_5)*M11+(dM11_5+dM33_5)*M22+(dM11_5+dM22_5)*M33
            -2*(dM12_5*M12+dM13_5*M13+dM23_5*M23))+(dM11_5*M22*M33+dM22_5*M11*M33+dM33_5*M11*M22)
            +2*(dM12_5*M13*M23+dM13_5*M12*M23+dM23_5*M12*M13)-(dM11_5*M23**2+dM22_5*M13**2+dM33_5*M12**2)
            -2*(dM12_5*M12*M33+dM13_5*M13*M22+dM23_5*M23*M11)

            dEp_6=mh2**2*(dM11_6+dM22_6+dM33_6)-mh2*((dM22_6+dM33_6)*M11+(dM11_6+dM33_6)*M22+(dM11_6+dM22_6)*M33
            -2*(dM12_6*M12+dM13_6*M13+dM23_6*M23))+(dM11_6*M22*M33+dM22_6*M11*M33+dM33_6*M11*M22)
            +2*(dM12_6*M13*M23+dM13_6*M12*M23+dM23_6*M12*M13)-(dM11_6*M23**2+dM22_6*M13**2+dM33_6*M12**2)
            -2*(dM12_6*M12*M33+dM13_6*M13*M22+dM23_6*M23*M11)

            dEp_7=mh2**2*(dM11_7+dM22_7+dM33_7)-mh2*((dM22_7+dM33_7)*M11+(dM11_7+dM33_7)*M22+(dM11_7+dM22_7)*M33
            -2*(dM12_7*M12+dM13_7*M13+dM23_7*M23))+(dM11_7*M22*M33+dM22_7*M11*M33+dM33_7*M11*M22)
            +2*(dM12_7*M13*M23+dM13_7*M12*M23+dM23_7*M12*M13)-(dM11_7*M23**2+dM22_7*M13**2+dM33_7*M12**2)
            -2*(dM12_7*M12*M33+dM13_7*M13*M22+dM23_7*M23*M11)

            dEp_8=mh2**2*(dM11_8+dM22_8+dM33_8)-mh2*((dM22_8+dM33_8)*M11+(dM11_8+dM33_8)*M22+(dM11_8+dM22_8)*M33
            -2*(dM12_8*M12+dM13_8*M13+dM23_8*M23))+(dM11_8*M22*M33+dM22_8*M11*M33+dM33_8*M11*M22)
            +2*(dM12_8*M13*M23+dM13_8*M12*M23+dM23_8*M12*M13)-(dM11_8*M23**2+dM22_8*M13**2+dM33_8*M12**2)
            -2*(dM12_8*M12*M33+dM13_8*M13*M22+dM23_8*M23*M11)
            
            ###dmhi=d(mh2)/d(pi).
            dmh1=-(dEz*dmz1+dEt*dtb1+dEm*dmu1+dEp_1)/dEh
            dmh2=-(dEz*dmz2+dEt*dtb2+dEm*dmu2+dEp_2)/dEh
            dmh3=-(dEz*dmz3+dEt*dtb3+dEm*dmu3+dEp_3)/dEh
            dmh4=-(dEz*dmz4+dEt*dtb4+dEm*dmu4+dEp_4)/dEh
            dmh5=-(dEz*dmz5+dEt*dtb5+dEm*dmu5+dEp_5)/dEh
            dmh6=-(dEz*dmz6+dEt*dtb6+dEm*dmu6+dEp_6)/dEh
            dmh7=-(dEz*dmz7+dEt*dtb7+dEm*dmu7+dEp_7)/dEh
            dmh8=-(dEz*dmz8+dEt*dtb8+dEm*dmu8+dEp_8)/dEh

            ###FT_mhi is the fine-tuning of mh2.
            FT_mh1=np.abs(mhu2/mh2*dmh1)
            FT_mh2=np.abs(mhd2/mh2*dmh2)
            FT_mh3=np.abs(ms2/mh2*dmh3)
            FT_mh4=np.abs(ak/mh2*dmh4)
            FT_mh5=np.abs(Yt/mh2*dmh5)
            FT_mh6=np.abs(al/mh2*dmh6)
            FT_mh7=np.abs(l/mh2*dmh7)
            FT_mh8=np.abs(k/mh2*dmh8)

            FT_mh={'FT_mh1':FT_mh1,'FT_mh2':FT_mh2,'FT_mh3':FT_mh3,'FT_mh4':FT_mh4,
            'FT_mh5':FT_mh5,'FT_mh6':FT_mh6,'FT_mh7':FT_mh7,'FT_mh8':FT_mh8}
            # print(FT_mh)
            return FT_mh
        FT_mh=mh_derivative()
        print(FT_mz,FT_mh)
        return FT_mz,FT_mh
    
    return mz_derivative() 
    

ZH_Finetuning(input_parameter)
