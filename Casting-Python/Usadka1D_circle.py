#Steps for a calculation:
#-Set parameters
#-CalcBcProp(Twat,Qwat)
#-CalcMatProp(C=0.1,....)
#-Set list for output
#-RunCalc()
from math import exp, tanh, log, pi
import sys
sys.path.append('C:\\Program Files\\ParaView 5.11.1\\bin\\Lib\\site-packages')
import vtk
class Usadka1D_circl:
    #---------- casting parameters---------------------------
    R=0.200/2     #radius of biller cross section, m    
    v=1.5           #casting speed, m/min
    dTemp=20       #temperature of overheating for liquid steel, K
    Tsr=15         #ambient temperature, Celsius
    MouldLevel=0.1 #mould level,m
    #---------- heat exchange parameters
    alfa_liq=2200  #htc for a border between solid and liquid steel, W/m2K
    lamda_liq=190  #conductivity in the area of luquid steel, W/m*K
    Zones=[[0.85,200.0],[0.925,100]] #(start of zone, m; water flux, l/m2*min)
    #---------- Output---------------------------
    Z_list=[0.1,0.15,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.801]
    Value_tst=[]     #list of test variables
    ValieName_tst=[] #list of names for test variables 
    #---------- mould parameters---------------------------
    Wat_thick=0.0035     #thickness of water layer, m
    Mould_lamda=370     #mould material conductivity, W/mK
    Mould_thick=0.012   #distance between water and mould surface, m
    Coat_lamda=80       #coating conductivity, W/mK
    Coat_thick0=0.0     #top thickness of coating, m
    Coat_thick1=0.0     #bottom thickness of coating, m
    Tflux_melt=1130     #melting temperature for mould flux, C
    alfa_flux_max=1700  #maximum HTC for flux, W/m2*K
    Rate_flux_drop=0.95  #coefficeint of HTC drop
    Rate_flux_melt=0.05 #coefficient of htan() defining temp range of melting
    Press_coeff=0.02   #coefficient of contact conductence
    #-----------steel properties----------------------
    lamda=29    #steel conductivity, W/m*K    
    L=272E+3    #heat of solidification, J/kg
    ro_sol=7300 #solid steel density, kg/m3
    ro_liq=7000 #liquid steel density, kg/m3
    Cl=795      #heat capacity for liquid steel, J/kg*K
    Cr=687      #heat capacity for solid steel, J/kg*K
    k=0.1413    #power for strain rate in the steel creep equation
    Tc=275      #thermal coefficient in the steel creep equation, K
    Ac=14546    #coefficient for stress in the steel creep equation, MPa
    #-----------variables for calculation
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    def CalcMatProp(self,C=0.0,Mn =0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S)                  #Solidus temperatere, Celsius
        self.Tlik=1536-(78*C+7.6*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S)  #Liquidus temperature, Celsius
        self.beta=(2.7-0.16*C+0.039*Mn-0.1*Si-0.019*Cr-0.016*Ni-0.5*P-0.25*S)/140000 #Thermal expansion coeff, 1/K
        self.Hl=(self.ro_liq*self.Cl+self.ro_sol*self.Cr)*(self.Tlik-self.Tsol)/2+(self.ro_liq+self.ro_sol)*self.L/2           #J/m3
        print('***Thermal Properties of steel:')
        print('Solidus temperature, Celsius: '+str(self.Tsol))
        print('Liquidus temperature, Celsius: '+str(self.Tlik))
        print('Expension coefficient *10^5: '+str(self.beta*1E+5))
        print(' ')
    def CalcBcProp(self,Twat,Qwat):
        #Twat - inlet water temeprature, Celsius
        #Qwat - water flow, l/min
        self.Qwat=Qwat
        self.Twat_in=Twat
        self.Twat=self.Twat_in     #current water temeprature
        lamda_wat=0.55748+0.0021525*Twat-0.0000097*Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*Twat)+2.52805091E-07 #m2/sec
        Pr=self.Prandtl(Twat)
        v_wat=Qwat/60000/pi/self.Wat_thick/(2*(self.R+self.Mould_thick)+self.Wat_thick)#m/sec
        Re=v_wat*2*self.Wat_thick/visc_wat
        self.alfa_wat0=0.023*lamda_wat/2/self.Wat_thick*Re**0.8*Pr**0.4*(1-0.45/(2.4+Pr))*(1+self.Wat_thick/(self.R+self.Mould_thick))**(0.16/Pr**0.15)
        print('***Heat transfer parameters:')
        print('Prandtl: '+str(Pr))
        print('Water speed, m/sec: '+str(v_wat))
        print('Reynolds: '+str(Re))
        print('Nominal HTC, kW/(m2K) :'+str(self.alfa_wat0/1000))
        print(' ')
    def HeatFlow(self,z,Ts):#W/m2
        if z<=self.Zones[0][0]:
            self.alfa_flux=self.v**0.8*(self.alfa_flux_max*(1-self.Rate_flux_drop/2*(1-tanh(self.Rate_flux_melt*(Ts-self.Tflux_melt))))+\
                self.Press_coeff*self.ro_liq*9.81*(z-self.MouldLevel))
            Coat_thick=self.Coat_thick0+(self.Coat_thick1-self.Coat_thick0)*z/self.Zones[0][0]
            self.alfa_wat=self.alfa_wat0*(1+self.Wat_thick/z)
            Q=(Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+Coat_thick)+\
                log(1+self.Mould_thick/(self.R+Coat_thick))/self.Mould_lamda+\
                log(1+Coat_thick/self.R)/self.Coat_lamda+\
                1/self.alfa_flux/self.R)
            TempW=self.Twat+Q*self.R/(self.R+self.Mould_thick)/self.alfa_wat
            self.alfa_wat*=(self.Prandtl(self.Twat_in)/self.Prandtl(TempW))**0.25
            return (Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+Coat_thick)+\
                    log(1+self.Mould_thick/(self.R+Coat_thick))/self.Mould_lamda+\
                    log(1+Coat_thick/self.R)/self.Coat_lamda+\
                    1/self.alfa_flux/self.R) #mould 
        else:
            iz=0
            while iz<len(self.Zones)-1 and self.Zones[iz+1][0]<z:iz+=1
            return 142/(5.5-z/30)*(Ts-self.Tsr)*self.Zones[iz][1]**0.55+5.670367E-8*((Ts+273)**4-(self.Tsr+273)**4) #spray
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value): #J/m3
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro_sol*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro_liq*self.Cl+self.Hl
        else: return (Value-self.Tsol)/(self.Tlik-self.Tsol)*self.Hl
    def Temperature(self,Value):
        if Value<0: return Value/self.ro_sol/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro_liq/self.Cl+self.Tlik
        else: return self.Tsol+(self.Tlik-self.Tsol)*Value/self.Hl
#-----------functions for shrinkage calculation
    def FindKsiC(self,j,i0,SpUs0,KsiZ0):
        ksi=[2*self.beta*self.SpT[j][i0]-KsiZ0/2-(KsiZ0*self.R*self.R/2+SpUs0*self.R)/i0/i0/self.dr/self.dr,\
            (KsiZ0*self.R*self.R/2+SpUs0*self.R)/i0/i0/self.dr/self.dr-self.beta*self.SpT[j][i0]-KsiZ0/2,\
            KsiZ0-self.beta*self.SpT[j][i0]]
        for i in range(i0,self.n):
            Value=3*self.beta*(self.SpT[j][i]+self.SpT[j][i+1])*(i+0.5)/i0/i0/2
            ksi[0]+=Value
            ksi[1]-=Value
        return ksi
    def NormalForce(self,j,n0,SpUs0,KsiZ0):
        Nf=[0,0]
        for i in range(n0,self.n+1):
            ksiC=self.FindKsiC(j,i,SpUs0,KsiZ0)
            KsiI=(2*((ksiC[0]-ksiC[1])**2+(ksiC[0]-ksiC[2])**2+(ksiC[1]-ksiC[2])**2))**0.5/3
            if KsiI==0: ValueLoc=0
            else:
                ValueLoc=2*self.Ac*KsiI**(self.k-1)*exp(-(self.T[j][i]+273)/self.Tc)/3
            if (i==n0)or(i==self.n):
                Nf[0]+=ValueLoc*(ksiC[1]-ksiC[0])*self.dr/2
                Nf[1]+=ValueLoc*(ksiC[2]-ksiC[0])*self.dr/2
            else:
                Nf[0]+=ValueLoc*(ksiC[1]-ksiC[0])*self.dr
                Nf[1]+=ValueLoc*(ksiC[2]-ksiC[0])*self.dr
        return Nf
    def SpeedUsadka(self,j,n0,SpUs0,KsiZ0):
        if SpUs0==0: dSpUs=-self.beta*100*self.R
        else: dSpUs=-100*self.Epsilon*SpUs0
        if KsiZ0==0: dKsiZ=-self.beta*100
        else: dKsiZ=-100*self.Epsilon*KsiZ0;
        eps=10*self.Epsilon
        while eps>self.Epsilon:
            F0=self.NormalForce(j,n0,SpUs0,KsiZ0)
            F1=self.NormalForce(j,n0,SpUs0+dSpUs,KsiZ0)
            F2=self.NormalForce(j,n0,SpUs0,KsiZ0+dKsiZ)
            if (F0[0]*F1[0])<0: dSpUs=dSpUs/2
            elif ((F1[0]-F0[0])*F0[0])>0: dSpUs=-dSpUs
            if (F0[1]*F2[1])<0: dKsiZ=dKsiZ/2
            elif ((F2[1]-F0[1])*F0[1])>0: dKsiZ=-dKsiZ
            eps=abs(2*dSpUs/(abs(SpUs0)+abs(dSpUs)))
            if (F0[0]*F1[0])<0: SpUs0+=dSpUs/2
            elif ((F1[0]-F0[0])*F0[0])>=0: SpUs0=SpUs0
            else: SpUs0+=dSpUs
            if (F0[1]*F2[1])<0: KsiZ0+=dKsiZ/2
            elif ((F2[1]-F0[1])*F0[1])>=0: KsiZ0=KsiZ0
            else: KsiZ0+=dKsiZ
        return SpUs0, KsiZ0
#---------------------------------------------------
    # n: number of divisions (node#0 - axis; node#n - outer surface of the billet)
    # kj: convergence coefficient for thermal calculations
    # Epsilon: Accuracy for calculations of shrinkage, if zero the calculation is turned off
    # CSVFile: Name of a file for output of general data
    # VTKFile: Name of vtk-file
    # Stiff: list of [Number of points, Max curvature rate, 1/(m*sec)]
    def RunCalc(self,n=100,kj=0.5,Epsilon=0.0, CSVFile='', VTKFile='', Stiff=None):        
        self.n=n
        self.dr=self.R/n
        dtau=kj*self.dr*self.dr*min(self.ro_liq*self.Cl,self.ro_sol*self.Cr)/4/max(self.lamda,self.lamda_liq)  #sek
        dZ=self.v*dtau/60   #m
        if CSVFile!='':
            f_csv=open(CSVFile,'w')
            f_csv.write('Axis, m; Shrinkage, mm; Surface Temperature, Celcius; Heat Flux, MVt/m2; Solid thickness,mm; ')
            f_csv.write('Liquid thickness, mm; Water temperature, Celsius; HTC water, kW/m2K; HTCflux, W/m2K\n')
        if Stiff!=None:
            f_stf=open('Stiffness.csv','w')
            f_stf.write('Z [m]; Stiffness [H*m2*sec] vs Curvature rate [1/(m*sec)]\n')
            for i in range(Stiff[0]): f_stf.write(';'+str(Stiff[1]*(i+1)/Stiff[0]))
            f_stf.write('\n')
        if VTKFile!='':
            mesh=vtk.vtkUnstructuredGrid()
            Points=vtk.vtkPoints()
            for Z in self.Z_list:
                for i in range(self.n+1):
                    Points.InsertNextPoint(i*self.dr,0,Z)
            mesh.Allocate((len(self.Z_list)-1)*self.n)
            mesh.SetPoints(Points)
            for j in range(len(self.Z_list)-1):
                for i in range(self.n):
                    mesh.InsertNextCell(vtk.VTK_QUAD,4,(j*(self.n+1)+i,j*(self.n+1)+i+1,(j+1)*(self.n+1)+i+1,(j+1)*(self.n+1)+i))
            vtkTemp=vtk.vtkFloatArray()
            vtkTemp.SetName('Temp')
            vtkTemp.SetNumberOfValues(len(self.Z_list)*(self.n+1))            
        print('_________________________________________')
        print(' Axis,m|Surf_T,C|Flux,kW/m2|Twat,C|Thickness,mm|')
        self.Epsilon=Epsilon
        self.SpT=[]
        self.T=[]
        self.T.append([])        
        out_i=0
        H=[]
        nk=n  # index of node after Tsol
        nF1=n # index of node after Tliq
        dRrate=0
        deltaR=0
        KsiZ=0
        DeltaF0=0
        DeltaF1=0
        Z=self.MouldLevel
        T0=self.Tlik+self.dTemp #initial temperature
        for i in range(0,n+1):
            self.T[0].append(T0)
            H.append(self.FuncTemp(T0)) #J/m3
        j=0
        while out_i<len(self.Z_list):
            Q=self.HeatFlow(Z,self.T[j][n])
    #--------------thickness---------------------
    #----- nk - index of node after Tsol
    #----- nF1 - index of node after Tliq
            if (nk<(n-2))and(nk>2):
                DeltaF0=((n-nk)*self.dr+(self.Tsol-self.T[j][nk])*self.dr/(self.T[j][nk-1]-self.T[j][nk]))
            if (nF1<(n-2))and(nF1>2):
                DeltaF1=((n-nF1)*self.dr+(self.Tlik-self.T[j][nF1])*self.dr/(self.T[j][nF1-1]-self.T[j][nF1]))
    #--------------Output------------------------------        
            if self.Z_list[out_i]-dZ/2<=Z<self.Z_list[out_i]+dZ/2:
                print(' {:6.3f}| {:6.1f} |  {:7.2f} | {:4.1f} |    {:5.1f}   |'.format(Z,self.T[j][n],Q/1000,self.Twat,DeltaF1*1000))
                if CSVFile!='':
                    #Axis, m; deltaR, mm; Surface Temperature, Celcius; Heat Flux, MVt/m2; Solid thickness,mm;
                    f_csv.write(str(Z)+';'+str(deltaR)+';'+str(self.T[j][n])+';'+str(Q/1000000)+';'+str(DeltaF0*1000)+';')
                    #Liquid thickness, mm; Water temperature, Celsius; HTC water, kW/m2K
                    f_csv.write(str(DeltaF1*1000)+';'+str(self.Twat)+';'+str(self.alfa_wat/1000)+';'+str(self.alfa_flux)+'\n')
                if Stiff!=None:
                    f_stf.write(str(Z))
                    for jj in range(Stiff[0]):
                        Jst=0
                        for i in range(nk,n+1):
                            if i==nk or i==n:
                                Jst+=(i**(2+self.k)*exp(-(self.T[j][i]+273)/self.Tc))*self.dr**(3+self.k)/2
                            else:
                                Jst+=(i**(2+self.k)*exp(-(self.T[j][i]+273)/self.Tc))*self.dr**(3+self.k)
                        f_stf.write(';'+str(Jst*self.Ac*(Stiff[1]*(jj+1)/Stiff[0])**(self.k-1)*(4-1.2091*self.k+0.4285*self.k**2)*10e+6))
                    f_stf.write('\n')
                if VTKFile!='':
                    for i in range(n+1):
                        vtkTemp.SetValue(out_i*(n+1)+i,self.T[j][i])
                out_i+=1
    #----------------------------------------------------
            if Z<=self.Zones[0][0]:self.Twat+=Q*dZ*2*pi*self.R/self.Qwat*60000/(4230000-3656.2*self.Twat-25.85*self.Twat**2)
            self.T.append([]) #j+1
            self.SpT.append([]) #j
    #----------Energy change calculation---------------
            for i in range(0,n):
#                if i<nF1-1:
#                    lamda=self.lamda_liq
#                elif i<nF1:
#                    lamda=self.alfa_liq*self.dr
#                else:
#                    lamda=self.lamda
                lamda=self.lamda
                if i==0: # axis
                    H[i]+=4*lamda*dtau/self.dr/self.dr*(self.T[j][i+1]-self.T[j][i])
                    H[i+1]+=-lamda*dtau/self.dr/self.dr*(1-1/2/(i+1))*(self.T[j][i+1]-self.T[j][i])
                else:
                    H[i]+=lamda*dtau/self.dr/self.dr*(1+1/2/i)*(self.T[j][i+1]-self.T[j][i])
                    if i==n-1: # outer surface
                        H[i+1]+=-lamda*dtau/self.dr/self.dr*(2-2/(4*n-1))*(self.T[j][i+1]-self.T[j][i])
                    else:
                        H[i+1]+=-lamda*dtau/self.dr/self.dr*(1-1/2/(i+1))*(self.T[j][i+1]-self.T[j][i])
            H[n]-=dtau*Q/self.dr*8*n/(4*n-1)
    #----------Temperature calculation for next level---------------
            for i in range(0,n+1):
                Value=self.Temperature(H[i])
                self.T[j+1].append(Value)
                self.SpT[j].append((Value-self.T[j][i])/dtau)
    #----------Shrinkage-----------------------------
            if Epsilon!=0:
                if self.T[j][n]<self.Tsol:
                    dRrate, KsiZ=self.SpeedUsadka(j,nk,self.beta*self.SpT[j][n]*self.R,KsiZ)
                deltaR+=dRrate*dtau*1000
    #------------------------------------------------
            nk=n
            nF1=n
            for i in range(0,n+1):
                if (nk==n)and(self.T[j+1][i]<self.Tsol): nk=i
                if (nF1==n)and(self.T[j+1][i]<self.Tlik): nF1=i
            Z+=dZ
            j+=1
        if CSVFile!='': f_csv.close()
        if Stiff!=None: f_stf.close()
        if VTKFile!='':
            mesh.GetPointData().SetScalars(vtkTemp)
            vtkoutput=vtk.vtkXMLUnstructuredGridWriter()
            vtkoutput.SetInputData(mesh)
            vtkoutput.SetFileName(VTKFile)
            vtkoutput.Write()
