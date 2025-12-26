def Log_message(text, file):
  print(txt)
  file.write(text+'\n')
#---------------------------------
class Solidification:
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
        self.ro_sol=7300    #solid steel density, kg/m3
        self.ro_liq=7000    #liquid steel density, kg/m3        
        self.Cl=680        # heat capacity for liquid steel, J/kg*K
        self.Cr=795        # heat capacity for solid steel, J/kg*K
        self.Hl=self.ro*((self.Cl+self.Cr)*(self.Tlik-self.Tsol)/2+self.L) #J/m3
        self.HTC1.ImportProp(self.lamda_liq, self.beta, self.kvis, self.thdif, self.Tsol,self.Cl,self.ro_liq)
        self.HTC2.ImportProp(self.lamda_liq, self.beta, self.kvis, self.thdif, self.Tsol,self.Cl,self.ro_liq) 
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value):
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro_sol*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro_liq*self.Cl+self.Hl
        else: return (Value-self.Tsol)/(self.Tlik-self.Tsol)*self.Hl
    def Temperature(self,Value):
        if Value<0: return Value/self.ro_sol/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro_liq/self.Cl+self.Tlik
        else: return self.Tsol+(self.Tlik-self.Tsol)*Value/self.Hl
