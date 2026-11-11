# Common style of the figures of the paper (serif fonts of at least 7 pt, TrueType fonts in PDF, 300 dpi)
import numpy as np, re, collections, json, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os; OUT=os.environ.get('FIGDIR','../results/figures/')   # where the manuscript figures go (default results/figures/)
os.makedirs(OUT,exist_ok=True); FIGDIR=OUT
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':0.6,'axes.edgecolor':'#555','axes.labelcolor':'#222','xtick.color':'#444','ytick.color':'#444','legend.frameon':False,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42})
COL={'BEOTTI':'#2a78d6','BEOUDS':'#eb6834','FRMAUB':'#1baf7a','BEHUMA':'#eda100','BEUCCE':'#e87ba4'}
MK={'BEOTTI':'o','BEOUDS':'s','FRMAUB':'^','BEHUMA':'D','BEUCCE':'v'}
NAME={'BEOTTI':'Ottignies','BEOUDS':'Oudsbergen','FRMAUB':'Maubeuge','BEHUMA':'Humain','BEUCCE':'Uccle'}
