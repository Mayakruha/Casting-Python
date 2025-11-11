from math import exp
#Steps for a calculation:
#-Create mould object
#-Set casting parameters
#-CalcMaterilaProp(C=0.0,....)
#-Set list for ouptut
#-RunCalc()
def Coating(z):# coating thickness, m
    return 0.0005+0.001*z/0.8
class Mould:
    def __init__(self,Height,Nch,Sch,Pch,Mould_thick,Thick,CoatFunc):
        #Nch - number of cooling channels
        #Sch - cross section of a cooling channel, m2
        #Pch - perimeter of a cooling channel, m
        self.Height=Height           # Mould height, m
        self.Twat_in=15              # Temperature of water [Celsius]
        self.Qwat=0                  # Water flow [l/min]
        self.Coat_lamda=80           # coating conductivity, W/mK
        self.Mould_lamda=370         # mould material conductivity, W/mK
        self.Taper=0.0005            # taper of a side over height, m
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Thb=Thick               # size for calculation of shrinkage
        self.deff=4*Sch/Pch          # effective diameter, m
        self.Wat_sec=Nch*Sch
        self.CoatFunc=CoatFunc
        self.flux_Tmelt=1115         # melting temperature for mould flux, C
        self.flux_Tliq=1145          # melting temperature for mould flux, C
        self.flux_lamda=1.5          # flux conductivity, W/mK
        self.flux_alfa_liq=5900      # HTC for luquid flux, W/m2*K
        self.flux_alfa_sol=1200      # HTC for solid flux, W/m2*K
        self.Zm=0
    def flux_alfa(self,T):
        if T>self.flux_Tliq: return self.flux_alfa_liq
        elif T<self.flux_Tmelt: return self.flux_alfa_sol
        else: return (T-self.flux_Tmelt)/(self.flux_Tliq-self.flux_Tmelt)*(self.flux_alfa_liq-self.flux_alfa_sol)+self.flux_alfa_sol
    def set_level(self,z):
        self.Coat_thick=self.CoatFunc(z)
        self.Rm=self.Mould_thick/self.Mould_lamda+self.Coat_thick/self.Coat_lamda
        self.z=z
        self.alfa_watz=self.alfa_wat0*(1+self.deff/z/2)
        if self.z>self.Zm:
            self.shrink=self.Thb*0.0085*(z-self.Zm)**0.5
    def mould_shape(self,z):#m
        return self.Taper*z/self.Height
    def HeatFlow(self,z,Ts):#W/m2
        if self.z==self.Zm:
            if Ts<=self.Tsol:
                self.flux_thickz=self.flux_thick_m
            else:
                self.flux_thickz=self.flux_thick_m*(1+2*(Ts-self.Tsol)/self.Tsol)
        elif self.z>self.Zm:
            self.flux_thickz=self.flux_thick_m+self.shrink-self.mould_shape(self.z)+self.mould_shape(self.Zm)
            if self.flux_thickz<0.0:self.flux_thickz=0.0
        Q=(Ts-self.Twat)/(1/self.alfa_watz+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(Ts)/self.v**0.8)
        self.TempW=self.Twat+Q/self.alfa_watz
        self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/(1/self.alfa_wat+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(Ts)/self.v**0.8) #mould 
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    # q - heat flux, W
    def HeatUp(self,q):
        self.Twat+=q/self.Qwat*60/(4230-3.6562*self.Twat-0.02585*self.Twat**2)
    def InPort(self): #should be replaced by OutPort from model
        return 0, 0, 0 #self.Zm, self.v, self.Tsol
    def Update(self):
        self.Zm, self.v, self.Tsol = self.InPort()
        print('\n** Heat transfer parameters')
        self.Twat=self.Twat_in                                       #current water temeprature
        lamda_wat=0.55748+0.0021525*self.Twat-0.0000097*self.Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*self.Twat)+2.52805091E-07 #m2/sec
        Pr=self.Prandtl(self.Twat)
        v_wat=self.Qwat/60000/self.Wat_sec #m/sec
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
class Usadka1D_slab:
    def __init__(self,mould):
        self.mould=mould
        self.mould.InPort=self.OutPort
        #---------- casting parameters---------------------------
        self.a=1.0/2        #half of slab width, m    
        self.b=0.155        #depth of a calculation (half of slab thickness), m
        self.v=0.95         #casting speed, m/min
        self.dTemp=0        #temperature of overheating for liquid steel, K
        self.Tsr=15         #ambient temperature, Celsius
        self.MouldLevel=0.1 #mould level,m
        self.alfa_liq=0     #HTC between liquid and solid steel (if b = half of slab thickness then alfa_liq = 0)
        self.Zones=[[mould.Height,200.0],[0.925,100]] #(start of zone, m; water flux, l/m2*min)
        #---------- Output---------------------------
        self.Z_list=[0.9,1.0,1.13,1.26,1.39,1.52]
        self.Zm=0
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn =0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S)                  #Solidus temperatere, Celsius
        self.Tlik=1536-(78*C+7.5*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S)  #Liquidus temperature, Celsius
        self.beta=(2.7-0.16*C+0.039*Mn-0.1*Si-0.019*Cr-0.016*Ni-0.5*P-0.25*S)/140000 #Thermal expansion coeff, 1/K
        self.lamda=29 #steel conductivity, W/m*K
        self.L=272E+3 #heat of solidification, J/kg
        self.ro=7200  #steel density, kg/m3 
        self.Cl=500   #heat capacity for liquid steel, J/kg*K
        self.Cr=680   #heat capacity for solid steel, J/kg*K
        self.k=0.1454 #power for strain rate in the steel creep equation
        self.Tc=266   #thermal coefficient in the steel creep equation, K
        self.Ac=11348 #coefficient for stress in the steel creep equation, MPa
        self.Hl=self.ro*((self.Cl+self.Cr)*(self.Tlik-self.Tsol)/2+self.L) #J/m3
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
    def NormalForce(self,j,n0,ksi0):
        Sum=0
        for i in range(n0,self.n+1):
            ValueLoc=ksi0-self.beta*self.SpT[j][i]
            if ValueLoc<0: ValueLoc=-((-ValueLoc) ** self.k)*exp(-(self.T[j][i]+273)/self.Tc)
            else: ValueLoc=(ValueLoc ** self.k)*exp(-(self.T[j][i]+273)/self.Tc)
            if (i==n0)or(i==self.n): Sum+=ValueLoc*self.dX/2
            else: Sum+=ValueLoc*self.dX
        return Sum
    def SpeedUsadka(self,j,n0,ksi0):
        if ksi0==0: dksi=-self.beta*100
        else: dksi=-100*self.Epsilon*ksi0
        eps=10*self.Epsilon
        while eps>self.Epsilon:
            ksi1=ksi0+dksi
            F0=self.NormalForce(j,n0,ksi0)
            F1=self.NormalForce(j,n0,ksi1)
            if (F0*F1)<0: dksi=dksi/2
            elif ((F1-F0)*F0)>0: dksi=-dksi
            elif (F0*F1)==0: dksi=0
            eps=abs(2*dksi/(abs(ksi0)+abs(ksi1)))
            if (F0*F1)<0: ksi0=(ksi1+ksi0)/2
            elif ((F1-F0)*F0)>=0: ksi0=ksi0
            else: ksi0=ksi1
        return ksi0
    def HeatFlow(self,z,Ts):#W/m2
        if z<=self.Zones[0][0]:
            return self.mould.HeatFlow(z,Ts) #mould 
        else:
            iz=0
            while iz<len(self.Zones)-1 and self.Zones[iz+1][0]<z:iz+=1
            return 142/(5.5-z/30)*(Ts-self.Tsr)*self.Zones[iz][1]**0.55+5.670367E-8*((Ts+273)**4-(self.Tsr+273)**4) #spray
    def OutPort(self):
        return self.MouldLevel, self.v, self.Tsol
    def RunCalc(self,n=300,kj=0.5,Epsilon=0.0001):
        # kj-convergence coefficient for thermal calculations
        self.n=n
        self.dX=self.b/n
        dtau=kj*self.dX*self.dX*self.ro*min(self.Cl,self.Cr)/4/self.lamda  #sek
        dZ=self.v*dtau/60   #m
        self.Epsilon=Epsilon
        self.SpT=[]
        self.T=[]
        self.T.append([])
        T0=self.Tlik+self.dTemp #initial temperature, K
        self.output=[]
        out_i=0
        H=[]
        H.append([])
        Q=[]
        nk=[]
        nk.append(n)
        ksi=0
        deltaB=0
#        SigmaX=0
        DeltaF0=0
        DeltaF1=0
        nF1=n
        Z=self.MouldLevel
        self.mould.Update()
        j=0
        for i in range(0,n+1):
            self.T[0].append(T0)
            H[0].append(self.FuncTemp(T0))
        while out_i<len(self.Z_list):
            if Z<=self.Zones[0][0]:
                self.mould.set_level(Z)
                if Z==self.mould.Zm:
                    Flag=False #for check of start of solidification
                    if self.T[j][n]>self.Tsol:Flag=True   
            Q.append(self.HeatFlow(Z,self.T[j][n]))
    #--------------thickness---------------------
            if (nk[j]<(n-2))and(nk[j]>2):
                DeltaF0=((n-nk[j])*self.dX+(self.Tsol-self.T[j][nk[j]])*self.dX/(self.T[j][nk[j]-1]-self.T[j][nk[j]]))
            if (nF1<(n-2))and(nF1>2):
                DeltaF1=((n-nF1)*self.dX+(self.Tlik-self.T[j][nF1])*self.dX/(self.T[j][nF1-1]-self.T[j][nF1]))
            Tliqavg=0
            for i in range(0,nF1):Tliqavg+=self.T[j][i]/nF1
    #--------------------------------------------        
            if self.Z_list[out_i]-dZ/2<=Z<self.Z_list[out_i]+dZ/2:
                #Axis, m; Shrinkage, mm; Surface Temperature, Celcius; Heat Flux, MVt/m2; Solid thickness,mm; Liquid thickness, mm; Temprature in the middle, Celsius; Water temp, Celsius 
                self.output.append((Z,deltaB,self.T[j][n],Q[j]/1000000,DeltaF0*1000,DeltaF1*1000,Tliqavg,self.mould.Twat))
                out_i+=1
            nk.append(n)
            nF1=n
            H.append([]) #j+1
            self.T.append([]) #j+1
            self.SpT.append([]) #j
            for i in range(0,n+1):
    #----------Temperature calculation---------------
                if i==0:
                    Value=2*self.lamda*dtau/self.dX/self.dX*(self.T[j][i+1]-self.T[j][i])+H[j][i]+2*dtau/self.dX*self.alfa_liq*(T0-self.T[j][i])
                elif i==n:
                    Value=2*self.lamda*dtau/self.dX/self.dX*(self.T[j][i-1]-self.T[j][i])+H[j][i]-2*dtau/self.dX*Q[j]
                else: Value=self.lamda*dtau/self.dX/self.dX*(self.T[j][i-1]+self.T[j][i+1]-2*self.T[j][i])+H[j][i]
                H[j+1].append(Value)
                Value=self.Temperature(Value)
                self.T[j+1].append(Value)
                if (nk[j+1]==n)and(Value<self.Tsol): nk[j+1]=i
                if (nF1==n)and(Value<self.Tlik): nF1=i
                self.SpT[j].append((Value-self.T[j][i])/dtau)
    #----------Shrinkage-----------------------------
            if self.T[j][n]<self.Tsol:
                ksi=self.SpeedUsadka(j,nk[j],self.beta*self.SpT[j][n])
#                SigmaX=self.Ac*(abs(2*(self.beta*self.SpT[j][n]-ksi)) ** self.k)*exp(-(self.T[j][n]+273)/self.Tc)
            deltaB+=ksi*self.a*dtau*1000
    #------------------------------------------------
            self.mould.HeatUp(Q[j]*dZ*2*self.a)
            Z+=dZ
            if Flag: self.mould.Zm=Z
            j+=1
