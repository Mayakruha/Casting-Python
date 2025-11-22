#Steps for a calculation:
#-Create objects for HTC
#-Create Casting1D object
#-CalcMaterilaProp(C=0.0,....)
#-Apply initial conditions
#-RunCalc()
#
import numpy as np
from math import exp
from meshio import Mesh, CellBlock
ScalarResList=['Time [sec]', 'minTemp [C]', 'Thickness [mm]', 'BulkTemp1 [C]', 'HTC1 [W]', 'Flux1 [W/m2]',
               'BulkTemp2 [C]', 'HTC2 [W]', 'Flux2 [W/m2]']
ArrayResList=['Temp [C]',]
#------------------------------------------
#----------Parent HTC class----------------
class HTC:  
    def ImportProp(self, lamda, beta, kvis, thdif, Tsol):
        self.lamda=lamda # conductivity of liqid steel, W/m*K
        self.beta=beta   # thermal expansion of luquid steel, 1/K
        self.kvis=kvis   # kinematic viscosity of luquid steel, m2/sec
        self.thdif=thdif # thermal diffusivity of luquid steel, m2/sec
        self.Tsol=Tsol   # Solidus temperatere, Celsius
    def htc(self, tm, temp, dtau): # tm - time, temp - temperature, dtau - time range
        return 0, 0
#
#----------Child HTC classes--------------
# HTC class to handle a table [time, htc]
class HTC_tab(HTC):
    def __init__(self, Tmet, tab):
        self.Tmet=Tmet   # function (tm) returns metal temparature over time
        self.htc1=tab    # table with two rows [time, htc]
    def htc(self, tm, temp, dtau):
        if self.Tmet(tm)<=temp:
            return 0, self.Tmet(tm)
        else:
            return np.interp(tm, self.htc1[0], self.htc1[1]), self.Tmet(tm)
# HTC class to handle convection
class HTC_conv(HTC):
    def __init__(self, Tmet, Length):
        self.Tmet=Tmet     # function (tm) returns metal temparature over time
        self.Length=Length # function (tm) returns length over time        
    def htc(self, tm, temp, dtau):
        if self.Tmet(tm)<=temp:
            return 0, self.Tmet(tm)
        else:
            Temp=self.Tmet(tm)
            L=self.Length(tm)
            Ra=9.81*self.beta*(Temp-temp)*L**3/self.kvis/self.thdif
            return 0.124*Ra**0.309*self.lamda/L, Temp
#
def Coating(z):# coating thickness, m
    return 0.0005+0.001*z/0.8
class CCM(HTC):
    def __init__(self, Height, Nch, Sch, Pch, Mould_thick, Thick, CoatFunc,
                 Mould_lamda=370, Coat_lamda=80):
        #Nch - number of cooling channels
        #Sch - cross section of a cooling channel, m2
        #Pch - perimeter of a cooling channel, m
        self.Height=Height           # Mould height, m
        self.deff=4*Sch/Pch          # effective diameter, m
        self.Wat_sec=Nch*Sch
        self.CoatFunc=CoatFunc
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Thb=Thick               # size for calculation of shrinkage
        self.Mould_lamda=Mould_lamda # mould material conductivity, W/mK
        self.Coat_lamda=Coat_lamda   # coating conductivity, W/mK
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    def flux_alfa(self,T):
        if T>self.flux_Tliq: return self.flux_alfa_liq
        elif T<self.flux_Tmelt: return self.flux_alfa_sol
        else: return (T-self.flux_Tmelt)/(self.flux_Tliq-self.flux_Tmelt)*(self.flux_alfa_liq-self.flux_alfa_sol)+self.flux_alfa_sol
    def SetParams(self, v, Zm, Qwat_mould, Taper, Qwat_spray, Twat_in=15, Tair=15,
                  flux_Tmelt=1115, flux_Tliq=1145, flux_lamda=1.5, flux_alfa_liq=5900, flux_alfa_sol=1200):
        self.v=v                          # Casting speed [m/min]
        self.Level=Zm                     # Level [m]
        self.Qwat_mould=Qwat_mould        # Water flow in the mould [l/min]
        self.Taper=Taper                  # taper of a side over height, m
        self.Twat_in=Twat_in              # Temperature of water [Celsius]
        self.Qwat_spray=Qwat_spray        # Water flow in the spray zone - list of (end of zone, m; water flux, l/m2*min)
        self.Tair=Tair                    # Temperature of air [Celsius]
        self.flux_Tmelt=flux_Tmelt        # melting temperature for mould flux, C
        self.flux_Tliq=flux_Tliq          # melting temperature for mould flux, C
        self.flux_lamda=flux_lamda        # flux conductivity, W/mK
        self.flux_alfa_liq=flux_alfa_liq  # HTC for luquid flux, W/m2*K
        self.flux_alfa_sol=flux_alfa_sol  # HTC for solid flux, W/m2*K
        self.Zm=Zm
        print('\n** Heat transfer parameters')
        self.Twat=self.Twat_in                                       #current water temeprature
        lamda_wat=0.55748+0.0021525*self.Twat-0.0000097*self.Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*self.Twat)+2.52805091E-07 #m2/sec
        Pr=self.Prandtl(self.Twat)
        v_wat=self.Qwat_mould/60000/self.Wat_sec #m/sec
        Re=v_wat*self.deff/visc_wat
        self.alfa_wat0=0.023*lamda_wat/self.deff*Re**0.8*Pr**0.4
        self.set_level(self.Zm)
        self.flux_thick_m=self.flux_lamda*((self.flux_Tmelt-self.Twat)/self.flux_alfa(self.Tsol)/self.v**0.8/(self.Tsol-self.flux_Tmelt)-self.Rm-1/self.alfa_watz)
        if self.flux_thick_m<0.0: self.flux_thick_m=0.0
        print('Prandtl: '+str(Pr))
        print('Water speed, m/sec: '+str(v_wat))
        print('Reynolds: '+str(Re))
        print('Nominal HTC, kW/(m2K): '+str(self.alfa_wat0/1000))
        print('Solid flux thickness at meniscus, mm: '+str(self.flux_thick_m*1000))
    def set_level(self, z):
        self.Coat_thick=self.CoatFunc(z)
        self.Rm=self.Mould_thick/self.Mould_lamda+self.Coat_thick/self.Coat_lamda
        self.alfa_watz=self.alfa_wat0*(1+self.deff/z/2)
    def mould_shape(self,z):#m
        return self.Taper*z/self.Height
    def htc(self, tm, temp, dtau):#W/m2
        z=self.Level+self.v*tm/60
        if z<=self.Height:
            if z>self.Zm:
                self.shrink=self.Thb*0.0085*(z-self.Zm)**0.5
            if temp>=self.Tsol:
                self.flux_thickz=self.flux_thick_m*(1+2*(temp-self.Tsol)/self.Tsol)
                self.Zm=z
            self.set_level(z)
            if z>self.Zm:
                self.flux_thickz=self.flux_thick_m+self.shrink-self.mould_shape(z)+self.mould_shape(self.Zm)
                if self.flux_thickz<0.0:self.flux_thickz=0.0
            Q=(temp-self.Twat)/(1/self.alfa_watz+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(temp)/self.v**0.8)
            self.TempW=self.Twat+Q/self.alfa_watz
            self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
            alfa=1/(1/self.alfa_wat+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(temp)/self.v**0.8)
            TempWat=self.Twat
            self.Twat+=dtau*(temp-self.Twat)*alfa/self.Qwat_mould*60/(4230-3.6562*self.Twat-0.02585*self.Twat**2)
            return alfa, TempWat #mould
        else: 
            iz=0
            while iz<len(self.Qwat_spray):
                if z>self.Qwat_spray[iz][0]:
                    iz+=1
                else:
                    break
            if iz>=len(self.Qwat_spray):
                alfa=5.670367E-8*((temp+273)**4-(self.Tair+273)**4)/(temp-self.Tair)
                return alfa, self.Tair
            else:
                alfa= 142/(5.5-z/30)*self.Qwat_spray[iz][1]**0.55+5.670367E-8*((temp+273)**4-(self.Tair+273)**4)/(temp-self.Tair) #spray
                return alfa, self.Tair
#--------------------------------------------------------
#---Functions applying initial conditions----------------
def meniscus(model, T0):
    for i in range(len(model.T)-1):
        model.T[i]=T0
    model.T[-1]=model.Tlik
def coldplate(model, Tcld, thick):
    k1=int((model.Dim-thick)/2/model.dX)-1
    k2=int((model.Dim+thick)/2/model.dX)+1
    alfa1, T1 = model.HTC1.htc(0, Tcld, 0)
    alfa2, T2 = model.HTC2.htc(0, Tcld, 0)
    for i in range(len(model.T)):
        if i<=k1:
            model.T[i]=T1
        elif i>=k2:
            model.T[i]=T2
        else:
            model.T[i]=Tcld
#------------------------------------------
#-------------Casting class----------------
class Casting1D:
    def __init__(self, HTC1, HTC2, Dim, Dim0=0, n=100, Radial=False):
        self.HTC1=HTC1       #htc class
        self.HTC2=HTC2       #htc class
        self.Dim=Dim
        self.n=n
        self.dX=(Dim-Dim0)/n
        self.T=np.zeros(n+1)
        self.RadFlg=Radial   # type of analysis. If True - radial
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn =0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0,N=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S)             #Solidus temperatere, Celsius
        self.Tlik=1535.1-(88*C+8*Si+5*Mn+1.5*Cr+4*Ni+5*Cu+3*Al+30*P+25*S+80*N)  #Liquidus temperature, Celsius
        self.beta=0.000031 # thermal expansion of luquid steel, 1/K
        self.kvis=7.91E-7  # kinematic viscosity of luquid steel, m2/sec
        self.thdif=5.45E-6 # thermal diffusivity of luquid steel, m2/sec
        self.lamda_liq=26  # conductivity of liqid steel, W/m*K
        self.lamda=32      # conductivity of solid steel, W/m*K
        self.L=272E+3      # heat of solidification, J/kg
        self.ro=7300       # steel density, kg/m3 
        self.Cl=680        # heat capacity for liquid steel, J/kg*K
        self.Cr=795        # heat capacity for solid steel, J/kg*K
        self.Hl=self.ro*((self.Cl+self.Cr)*(self.Tlik-self.Tsol)/2+self.L) #J/m3
        self.HTC1.ImportProp(self.lamda_liq, self.beta, self.kvis, self.thdif, self.Tsol)
        self.HTC2.ImportProp(self.lamda_liq, self.beta, self.kvis, self.thdif, self.Tsol)        
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value):
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro*self.Cl+self.Hl
        else: return self.ro*(Value-self.Tsol)/(self.Tlik-self.Tsol)*(self.L+(self.Tlik-self.Tsol)*self.Cr+(self.Cl-self.Cr)*(Value-self.Tsol)/2)
    def Temperature(self,Value):
        if Value<0: return Value/self.ro/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro/self.Cl+self.Tlik
        else:
            Tk=Value/self.Hl*(self.Tlik-self.Tsol)+self.Tsol
            eps=10*self.Epsilon
            while eps>self.Epsilon:
                dH=Value-self.FuncTemp(Tk)
                eps=abs(dH/self.Hl)
                Tk=dH/self.ro/self.Cef(Tk)+Tk
            return Tk
    def RunCalc(self, FullTime, kj=0.5, Epsilon=0.0001, out_dtau=1):
        # kj-convergence coefficient for thermal calculations
        self.dtau=kj*self.dX*self.dX*self.ro*min(self.Cl,self.Cr)/4/self.lamda  #sek
        self.Epsilon=Epsilon
        self.results={}
        for Name in ScalarResList+ArrayResList:
            self.results[Name]=[]
        out_time=0
        iter_time=0
        minTemp=min(self.T)
        n=len(self.T)
        H=np.zeros(n)
        for i in range(n):
            H[i]=self.FuncTemp(self.T[i])
        print('Solidus Temperature: {:6.1f}'.format(self.Tsol))
        print('Liquidus Temperature: {:6.1f}'.format(self.Tlik))
        print('********************************************')
        print(' Time, sec | Min Temp, C | Thickness, mm | Bulk Temp, C | HTC 1, kW/m2K | HTC 2, kW/m2K')
        while iter_time<FullTime and minTemp<self.Tlik:
    #--------------------------------------------
            k1=0
            k2=0
            for i in range(n):
                if (k1==0)and(self.T[i]<self.Tlik):
                    k1=i-1
                if self.T[i]<=self.Tlik:
                    k2=i+1                
            alfa1, T1 = self.HTC1.htc(iter_time,self.T[k1], self.dtau)
            alfa2, T2 = self.HTC2.htc(iter_time,self.T[k2], self.dtau)
            Q1=alfa1*(T1-self.T[k1])
            Q2=alfa2*(T2-self.T[k2])
            if iter_time>=out_time:
                self.results['Time [sec]'].append(iter_time)
                self.results['Temp [C]'].append(self.T.copy())
                self.results['minTemp [C]'].append(minTemp)
                self.results['Thickness [mm]'].append(self.dX*(k2-k1-2)*1000)
                self.results['BulkTemp1 [C]'].append(T1)
                self.results['HTC1 [W]'].append(alfa1)
                self.results['Flux1 [W/m2]'].append(Q1)
                self.results['BulkTemp2 [C]'].append(T2)
                self.results['HTC2 [W]'].append(alfa2)
                self.results['Flux2 [W/m2]'].append(Q2)
                print('  {:6.1f}   |    {:6.1f}   |    {:6.2f}     |    {:6.1f}    |    {:6.3f}     |    {:6.3f}'.format(iter_time, minTemp, self.dX*(k2-k1-2)*1000, T1, alfa1/1000, alfa2/1000))
                out_time+=out_dtau
    #----------Temperature calculation---------------
            H1=self.FuncTemp(T1)
            H2=self.FuncTemp(T2)
            for i in range(n):
                if i<k1:
                    H[i]=H1
                elif i>k2:
                    H[i]=H2
                elif i==k1:
                    H[i]=2*self.lamda*self.dtau/self.dX/self.dX*(self.T[i+1]-self.T[i])+H[i]+2*self.dtau/self.dX*Q1
                elif i==k2:
                    H[i]=2*self.lamda*self.dtau/self.dX/self.dX*(self.T[i-1]-self.T[i])+H[i]+2*self.dtau/self.dX*Q2
                else:
                    H[i]=self.lamda*self.dtau/self.dX/self.dX*(self.T[i-1]+self.T[i+1]-2*self.T[i])+H[i]
            for i in range(n):
                if i<k1:
                    self.T[i]=T1
                elif i>k2:
                    self.T[i]=T2
                else:
                    self.T[i]=self.Temperature(H[i])
            iter_time+=self.dtau
            minTemp=min(self.T)
#---------------------------------
#--------FUNCTIONS----------------
#---------------------------------
def output_csv(FileName, model):
    f=open(FileName,'w')
    for Name in ScalarResList:
        f.write(Name+';')
    f.write('\n')
    n=len(model.results[Name])
    for i in range(n):
        for Name in ScalarResList:
            f.write(str(model.results[Name][i])+';')
        f.write('\n')
    f.close()
#---------------------------------
def output_solid_vtu(FileName, model, Z):
    Dict_nodes={}
    Dict_nodes[0]={}
    coords=[]
    temp=[]
    cells=[]
    node=0
    for i in range(len(model.results['Temp [C]'])-1):
        Dict_nodes[(i+1)]={}
        for j in range(len(model.results['Temp [C]'][i])-1):
            if model.results['Temp [C]'][i][j]<model.Tlik and model.results['Temp [C]'][i][j+1]<model.Tlik and model.results['Temp [C]'][(i+1)][j]<model.Tlik and model.results['Temp [C]'][(i+1)][j+1]<model.Tlik:
                for i1 in range(2):
                    for j1 in range(2):
                        if not j+j1 in Dict_nodes[(i+i1)]:
                            Dict_nodes[(i+i1)][j+j1]=node
                            coords.append(((j+j1)*model.dX*1000,Z(model.results['Time [sec]'][i+i1])*1000,0))
                            temp.append(model.results['Temp [C]'][(i+i1)][j+j1])
                            node+=1
                cells.append((Dict_nodes[i][j],Dict_nodes[(i+1)][j],Dict_nodes[(i+1)][j+1],Dict_nodes[i][j+1]))
    model_out=Mesh(coords,[CellBlock('quad',cells),],point_data={'Temp [C]':temp})
    model_out.write(FileName)
