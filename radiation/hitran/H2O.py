
from hapi import *

db_begin('data')

fetch('CO2',2,1,0,10000)

for i in range (0,52):
	a = i*10 + 500
	b = "specCO2/%dK.txt" % (a)	
	print "doing temperature: %d" % (a)
	
	nu,coef = absorptionCoefficient_Lorentz(SourceTables='CO2',Environment={'T':a,'p':1},GammaL='gamma_self',HITRAN_units=False,File=b)
