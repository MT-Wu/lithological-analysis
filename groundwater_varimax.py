# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
'import matplotlib.pyplot as plt'

from starpy.stats.eof import srpca

data = np.genfromtxt('waterlevel.csv', dtype=float, delimiter=',')

dataT = np.transpose(data, axes=None)

m = 5  
"""EOF1~5"""
NewEOFs,NewEC,Newlambda=srpca(dataT,20,1)
""" 解釋比例加到20就是總和100% """

SUMlambda=np.sum(Newlambda)
percent=Newlambda/SUMlambda*100

MEOFs=NewEOFs[:,0:m]
MECs=NewEC[:,0:m]
Mpercent=percent[0:m]

col=[]
col2=[]
for i in xrange(1,m+1):
  col.append('EOF%d' % i)
  col2.append('EC%d' % i)
  

EOFdf=pd.DataFrame(MEOFs,columns=col)
ECdf=pd.DataFrame(MECs,columns=col2)
"percentdf=pd.DataFrame(Mpercent,columns=col)"



EOFdf.to_csv('./ROT_Test_EOF_temporal_norm.csv',index='FALSE')
ECdf.to_csv('./ROT_Test_EC_spatial_norm.csv',index='FALSE')
"Mpercent.pd.to_csv('./ROT_Taichung_percent_20151002.csv',index='FALSE')"


