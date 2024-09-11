#!/usr/bin/env python
import sys
from math import fabs, isnan, log,sqrt
from scipy import stats
import numpy as np
from numpy import log2, log10
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import signal
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

its = "Intensity "
lfq = "LFQ intensity "
q_dir = sys.argv[4]
prefix = sys.argv[5:]

scale = 1.5
bg_noise_cut = 5.0
bg_tol = 1e-2

#load fasta to get kD data
lib_fn = sys.argv[1]
map_mass_kD = {}
for record in SeqIO.parse(lib_fn, "fasta"):
  try:
    analysed_seq = ProteinAnalysis(str(record.seq))
    pro_name = record.id.split('|')[1]
    map_mass_kD[pro_name] = int(analysed_seq.molecular_weight()/1000)
    #print(pro_name, map_mass_kD[pro_name])
  except:
    #print "Warning:", pro_name, "has Z"
    pass

#markers
markers = open( sys.argv[3], 'r' ).readlines()
ms = []
ls = []
for marker in markers:
  es = marker.split()
  m = float(es[0])
  l = float(es[1])
  ms.append(m)
  ls.append(l)
ms = np.array(ms)
ls = np.array(ls)

print("#marker:", ls, ms)

def decay( x, p ):
  return p[0] / ( x + p[1] ) ** scale + p[2]

def residuals( p, y, x ):
  return (y - decay(x, p))**2

#fit N_frac, kD function
p0 = [ 50000.0, 50.0, -10.0 ]
plsq = leastsq(residuals, p0, args=(ls, ms))
a = plsq[0][0]
b = plsq[0][1]
c = plsq[0][2]
print( "DB:", a, b, c )

pdf = matplotlib.backends.backend_pdf.PdfPages("mass_label.pdf")
f, ax = plt.subplots(1, 1)
ax.scatter(ls, ms, marker="o", s=40, edgecolors='b')
fit_l = range( 1, 100 )
fit_mass = [  a / (b+l) ** scale + c for l in fit_l ]
ax.plot(fit_l, fit_mass)
ax.set_xlabel('Cutting marker (mm)')
ax.set_ylabel('Protein mass (kD)')
pdf.savefig(f)
pdf.close()

#fraction list
experiment_lst = []
map_frac_mass = {} #start from 1
frac_xtics = []
lines = open(sys.argv[2], 'r').readlines()
for n, l in enumerate(lines):
  es = l.strip().split()
  experiment_lst.append(es[0]) 
  #label
  l1 = float( es[1] )
  l2 = float( es[2] )
  #mass
  m1 = a / (b+l1) ** scale + c
  m2 = a / (b+l2) ** scale + c
  if l1<=0 or m1>1000.0: m1 = 1000.0
  map_frac_mass[n+1] = ( m2, m1 )
  frac_xtics.append( str(int(m2)) + '-' + str(int(m1)) )

#fit multiple gaussian?
def gau( x, *params ):
  N = int(len( params ) / 3)
  sgau = np.zeros(x.size)
  for i in range(N):
    a = params[i*3]
    x0 = params[i*3+1]
    sigma = params[i*3+2]
    sgau += a*np.exp( -(x-x0)**2/(2*sigma**2) )
  return sgau

def mask_peak( I, ndx ):
  N = len(I)
  I_msk = [0.0]*N
  I_msk[ndx] = I[ndx]
  i = ndx+1
  while i<N-1 and i<=ndx+2:
    I_msk[i] = I[i]
    if I[i]<I[i-1] and I[i]<I[i+1]: break
    i += 1
  i = ndx-1
  while i>0 and i>=ndx-2:
    I_msk[i] = I[i]
    if I[i]<I[i-1] and I[i]<I[i+1]: break
    i -= 1
  #print("DB mask:", I_msk)
  return I_msk

def analyze_gussian_peaks( gau_params ):
  #for reference, output mono rate
  #for other sample, output mono rate, largest shift
  if len(gau_params) == 0: return
  if len(gau_params[0]) == 0: return
  if len(gau_params[0][0]) == 0: return
  x0 = gau_params[0][0][0]
  #print("mono peak x0", x0)
  #print(gau_params)
  for peaks in gau_params:
    mono_prob = 0.0
    sum_prob = 0.0
    min_x = 99 #small one
    mono_x = x0 #big one
    mono_prob = 0
    for x, A, s in peaks:
      if x < x0 + 1.0:
        #skip degration peak
        area = A * s
        sum_prob += area
        if mono_prob == 0:
          mono_prob = area
          mono_x = x

        #skip small peak
        if x < min_x and area>1e-6:
          min_x = x
    if sum_prob>0:
      print("rate", mono_prob/sum_prob, "shift", mono_x - min_x, end=" ")
  return

lines = open(q_dir+"/proteinGroups.txt", 'r').readlines()
#tags
tags = {}
es = lines[0].strip().split('\t')
for i, e in enumerate(es):
  tags[e] = i

def process( name, *Is ):
  print("PRO:", name)
  Nsample = len(Is[0])
  Nfrac = len(Is[0][0])
  print(Nsample, Nfrac)
  mass_lst = []
  for p in name.split(';'):
    if p in map_mass_kD.keys():
      mass_lst.append( map_mass_kD[p] )
  if len(mass_lst) == 0: 
    print(p, "No mass found!")
    return
  # mono mass of the protein (group)
  guess_mass = np.median(mass_lst)
  print("guess_mass", guess_mass)
  # the closest marker of that mass 
  guess_mark = 0
  for i in range(1, Nfrac+1):
    m1, m2 = map_frac_mass[i]
    if m2>=guess_mass and m1<=guess_mass: #m1<r<m2
      guess_mark = i
      break
  if guess_mark == 0:
    if guess_mass < 20: guess_mark=Nfrac
    if guess_mass > 100: guess_mark=1
  print("guess_mark:", guess_mark, guess_mass)

  mono_peak_ndx = np.argmax(Is[0][0]) 
  mono_peak_save = -1
  mono_maxsig_save = -1

  print("mono_frac:", mono_peak_ndx+1)
  if mono_peak_ndx+1 == 1 and guess_mark == 1: return
  elif mono_peak_ndx+1 == Nfrac and guess_mark == Nfrac: return

  tnames = [ t for t in name.split(';') if "CON_" not in t ]
  xtics = np.arange(Nfrac)+1

  print("Nsample=", Nsample)
  f, ax = plt.subplots(1, Nsample, figsize=(4*Nsample, 4.8))
  plt.gca().invert_yaxis()
  #mask guess peak
  ax[0].axhline( guess_mark+1, -2, 4, linestyle='-', linewidth=1, c='black')
  ax[0].set_title("\n".join(tnames))

  gau_params = []
  for i, I in enumerate(Is[0]):
    #adjust
    save_mono_peak_ndx = mono_peak_ndx
    shift = []
    true_shift = 0
    if mono_peak_ndx>0: shift.append(-1)
    if mono_peak_ndx<Nfrac-1: shift.append(1)
    max_signal = I[mono_peak_ndx]
    for j in shift:
      if I[mono_peak_ndx+j]>max_signal:
        max_signal = I[mono_peak_ndx+j]
        true_shift = j
    mono_peak_ndx += true_shift

    #for each sample
    ax[i].axvline(0,0,Nfrac,linestyle='--',linewidth=1,c='black')
    for j in range(Nfrac):
      ax[i].axhline(j+1, -5, 5, linestyle='--', linewidth=0.4, c='grey', alpha=0.5)
    ax[i].set_ylim([Nfrac, 0])
    ax[i].barh( xtics-0.5, I, 1.0, color='red', alpha=0.32 )

    #fit gaussian
    final_params = []
    final_peaks = []
    #first peak
    p0 = [ max_signal, mono_peak_ndx+1.0, 1.0 ] 
    b0 = [(max_signal*0.9, mono_peak_ndx+0.5, 0.2), (max_signal/0.9, mono_peak_ndx+1.5, 4.0)]
    try:
      if mono_peak_save == -1:
        #first
        popt, pcov = curve_fit(gau, xtics, I, p0, bounds=b0)
        #popt, pcov = curve_fit(gau, xtics, mask_peak(I, mono_peak_ndx), p0, bounds=b0)
      else:
        popt, pcov = curve_fit(gau, xtics, mask_peak(I, mono_peak_ndx), p0, bounds=b0)
        #popt, pcov = curve_fit(gau, xtics, I, p0, bounds=b0)
      top = popt[0]
      mu0 = popt[1]
      sig = popt[2]
      final_params = list(popt)
    except:
      print("Fitting failed (1st)")
      #continue
      top = 1.0
      mu0 = mono_peak_ndx + 1.0
      sig = 1.0
      final_params = [mu0, top, sig]
    final_peaks.append( (mu0, top, sig) )
    ax[i].text(top, mu0+0.1, '*', fontweight='bold', fontsize=14 )
    #print("DB: (", mono_peak_ndx, mono_peak_ndx+2, ")")
    if mono_peak_save == -1:
      #first
      mono_peak_save = mu0
      mono_maxsig_save = top
    print("fitting 1:", final_params)

    #second peak
    I2 = []
    for x, y in zip( xtics, I ):
      if abs(x-mu0)<=1.6:
        I2.append(bg_tol)
      else:
        p_all = gau( x, *final_params )
        y -= p_all[0]
        if y<bg_tol: y=bg_tol
        I2.append(y)
    #print("DB I2:", I2)
    max_peak_ndx = np.argmax(I2) 
    #print("DB: max mu0", max_peak_ndx, mu0)
    max_signal = I2[max_peak_ndx]
    bg_signal = gau( max_peak_ndx, *final_params )
    if bg_signal<bg_tol: bg_signal = bg_tol
    if max_signal<bg_tol: max_signal = bg_tol
    try:
      sn = log(max_signal)-log(bg_signal)
    except:
      sn = 0
    #print( "DB: sn=", sn )
    #print( "DB: max_signal=", max_signal, "save=", mono_maxsig_save )
    if sn > bg_noise_cut and max_signal>mono_maxsig_save/50:
      p0 = [ max_signal, max_peak_ndx+1.0, 1.0 ]
      b0 = [( max_signal*0.9, max_peak_ndx+0.5, 0.2 ), ( max_signal/0.9, max_peak_ndx+1.5, 4.0 ) ]
      try:
        #popt, pcov = curve_fit(gau, xtics, I2, p0, bounds=b0)
        popt, pcov = curve_fit(gau, xtics, mask_peak(I2, max_peak_ndx), p0, bounds=b0)
      except:
        print("Fitting failed (2nd)")
        continue
      final_params += list(popt)
      final_peaks.append( (popt[1], popt[0], popt[2]) ) 
      top = popt[0]
      mu0 = popt[1]
      ax[i].text(top, mu0+0.1, '**', fontweight='bold', fontsize=12 )
      print("fitting 2:", final_params)
    else:
      #draw curve
      x_fit = np.linspace( 0, Nfrac, 100 )
      y_fit = gau( x_fit, *final_params )
      ax[i].plot( y_fit, x_fit-0.5, color='blue' )
      gau_params.append( final_peaks )
      continue

    #third peak
    I3 = []
    for x, y in zip( xtics, I ):
      if abs(x-1-mono_peak_ndx)<=1 and abs(x-max_peak_ndx)<=1:
        I3.append(bg_tol)
      else:
        p_all = gau( x, *final_params )
        y -= p_all[0]
        if y<bg_tol: y=bg_tol
        I3.append(y)
    #print("DB I3:", I3)
    max_peak_ndx = np.argmax(I3) 
    max_signal = I3[max_peak_ndx]
    bg_signal = gau( max_peak_ndx, *final_params )
    if bg_signal<bg_tol: bg_signal=bg_tol
    if max_signal<bg_tol: max_signal=bg_tol
    try:
      sn = log(max_signal)-log(bg_signal)
    except:
      sn = 0
    if sn > bg_noise_cut and max_signal>mono_maxsig_save/50:
      p0 = [ max_signal, max_peak_ndx+1.0, 1.0 ]
      b0 = [( max_signal*0.9, max_peak_ndx+0.5, 0.2 ), ( max_signal/0.9, max_peak_ndx+1.5, 4.0 ) ]
      try:
        popt, pcov = curve_fit(gau, xtics, I3, p0, bounds=b0)
      except:
        print("Fitting failed (3rd)")
        continue
      final_params += list(popt)
      final_peaks.append( (popt[1], popt[0], popt[2]) ) 
      top = popt[0]
      mu0 = popt[1]
      ax[i].text(top, mu0+0.1, '***', fontweight='bold', fontsize=12 )
      print("fitting 3:", final_params)
    else:
      #draw curve
      x_fit = np.linspace( 0, Nfrac, 100 )
      y_fit = gau( x_fit, *final_params )
      ax[i].plot( y_fit, x_fit-0.5, color='blue' )
      gau_params.append( final_peaks )
      continue

    #4 peak
    I4 = []
    for x, y in zip( xtics, I ):
      if abs(x-1-mono_peak_ndx)<=1 and abs(x-max_peak_ndx)<=1:
        I4.append(bg_tol)
      else:
        p_all = gau( x, *final_params )
        y -= p_all[0]
        if y<bg_tol: y=bg_tol
        I4.append(y)
    max_peak_ndx = np.argmax(I4) 
    max_signal = I4[max_peak_ndx]
    bg_signal = gau( max_peak_ndx, *final_params )
    #print("DB I4:", I4, max_peak_ndx, bg_signal)
    if bg_signal<bg_tol: bg_signal=bg_tol
    if max_signal<bg_tol: max_signal=bg_tol
    try:
      sn = log(max_signal)-log(bg_signal)
    except:
      sn = 0
    if sn > bg_noise_cut and max_signal>mono_maxsig_save/50:
      p0 = [ max_signal, max_peak_ndx+1.0, 1.0 ]
      b0 = [( max_signal*0.9, max_peak_ndx+0.5, 0.2 ), ( max_signal/0.9, max_peak_ndx+1.5, 4.0 ) ]
      try:
        popt, pcov = curve_fit(gau, xtics, I4, p0, bounds=b0)
      except:
        print("Fitting failed (4th)")
        continue
      final_params += list(popt)
      final_peaks.append( (popt[1], popt[0], popt[2]) ) 
      top = popt[0]
      mu0 = popt[1]
      ax[i].text(top, mu0+0.1, '****', fontweight='bold', fontsize=12 )
      print("fitting 4:", final_params)
    else:
      #draw curve
      x_fit = np.linspace( 0, Nfrac, 100 )
      y_fit = gau( x_fit, *final_params )
      ax[i].plot( y_fit, x_fit-0.5, color='blue' )
      gau_params.append( final_peaks )
      continue

    #draw curve
    x_fit = np.linspace( 0, Nfrac, 100 )
    y_fit = gau( x_fit, *final_params )
    ax[i].plot( y_fit, x_fit-0.5, color='blue' )
    gau_params.append( final_peaks )

    mono_peak_ndx = save_mono_peak_ndx

  ##########################################
  #statistics
  ##########################################
  print("STAT", name, end=" ")
  analyze_gussian_peaks( gau_params )
  print()
  ##########################################

  #finish
  pdf.savefig(f)
  plt.close('all')
  return

pdf = matplotlib.backends.backend_pdf.PdfPages("compare-outputs.pdf")

for l in lines[1:]:
  es = l.strip().split('\t')
  pros = es[tags["Protein IDs"]]
  
  #load 
  samples = []
  sum_lfq = 0.0
  for p in prefix:
    print("SAMPLE:", p)
    lfq_intens = []
    for exp in experiment_lst:
      lfq_intens.append( float(es[tags[lfq+p+exp]]) )
      sum_lfq += lfq_intens[-1]
    print(lfq_intens)
    samples.append(lfq_intens)
  if sum_lfq > 10.0:
    process(pros, samples)

pdf.close()

