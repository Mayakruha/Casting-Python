from FEMtoolkit import import_abq
model=import_abq('Sec_410_orig.inp')
model.NodeIntoSurf('Outer')
mesh=model.ReNumb()
mesh.NSets={}
mesh.export_abq('D://VTP//CCM//3DCCM//Sec_410.inp')
