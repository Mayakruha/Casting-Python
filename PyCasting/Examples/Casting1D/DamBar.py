from Casting1D import Casting1D, HTC_tab, HTC_conv, output_csv, output_solid_vtu, coldplate
import numpy as np
#
tm_full=830
Tinit=1570  #Tinit=1570
thick=0.014 # thickness, m
#
def Tmet(tm):
    return Tinit-tm*0.0482
def Depth(tm):
    if tm<25: return 15-0.4
    elif tm<75: return 14.6-(tm-25)*0.3/60/50
    elif tm<125: return 14.35-(tm-75)*0.7/60/50
    elif tm<245: return 13.767-(tm-125)/120
    elif tm<830: return 11.789-(tm-245)*1.2/60
    else: return 0
#-----------HTC----------------------
tab1=np.array([[  0.0,36.9,40.9,44.9,50.9,60.9,76.9,82.9,100.9,111.9,131.9,150.9,169.9,201.4,210.4,220.4,225.4,259.9,830],
                [12054,5059,5015,4845,4252,3940,3140,3726, 1938, 4068, 4995, 4507, 3702, 1813, 2398, 2932, 3007, 2912,  0]])
#--------------------------------------------
htc1=HTC_tab(Tmet, tab1)
htc2=HTC_conv(Tmet, Depth)
model=Casting1D(htc1, htc2, thick*1.8, n=150)
model.CalcMaterialProp(C=0.22, Mn =0.65, Si=0.3, P=0.04, S=0.05, Cu=0.3, Ni=0.3, Cr=0.3, N=0.008)
coldplate(model, 30, thick)
model.RunCalc(tm_full)
#
output_csv('Result_14.csv', model)
output_solid_vtu('Result_14.vtu', model, lambda tm: tm)
