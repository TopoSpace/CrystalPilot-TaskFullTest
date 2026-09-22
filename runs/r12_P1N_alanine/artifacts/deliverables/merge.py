"""Equivalent-reflection scaling, orthorhombic merging and extinction tests."""
import json,itertools
import numpy as np
from scipy.optimize import least_squares
from raw_frames import OUT

def keys_h(h,kind='mmm'):
 if kind=='mmm':return np.abs(h).astype(int)
 if kind=='-1':
  h=h.astype(int).copy()
  for i,row in enumerate(h):
   nz=row[row!=0]
   if len(nz) and nz[0]<0:h[i]*=-1
  return h

def merge():
 a=np.concatenate([np.load(OUT/f'integrated_run{r}.npy') for r in range(1,7)])
 valid=(a[:,10]==0)&(a[:,17]>=2)&(a[:,8]>.025)&np.isfinite(a[:,11])&(a[:,12]>0)
 masks=np.load(OUT/'static_masks.npy')
 shadow=masks[a[:,0].astype(int)-1,np.rint(a[:,6]).astype(int),np.rint(a[:,5]).astype(int)]
 np.savetxt(OUT/'masked_observations.csv',a[shadow],delimiter=',',header='run,h,k,l,omega,x,y,module,Linv,P,edge,I,sigma,Iraw,varraw,counts,bgsum,nframes,wom,wpos,wom2',comments='',fmt='%.8g')
 print('Raw-background shadow-mask exclusions',int(np.sum(valid&shadow)),flush=True)
 valid&=~shadow
 data=a[valid];h=keys_h(data[:,1:4]);uniq,inv=np.unique(h,axis=0,return_inverse=True);ng=len(uniq)
 raw=data[:,11];sig=data[:,12]
 # Expected Poisson variance: background contribution plus the group's signal,
 # not each observation's own signal. Prevent low counts from getting high weight.
 bvar=np.maximum(data[:,14]-data[:,13],0)
 lc=data[:,8]/data[:,9]
 def expected_sigma(sc,groupmean,inv):
  c=lc*sc
  return np.sqrt(np.maximum(bvar*c*c+np.maximum(groupmean[inv],0)*c,1e-12))
 # Only symmetry equivalents determine scale, never a structural model.
 X=np.column_stack([(data[:,0]==r).astype(float) for r in range(2,7)]+[(data[:,7]==1).astype(float)])
 groupids=[np.where(inv==i)[0] for i in range(ng)]
 count0=np.bincount(inv,minlength=ng)
 pscale=np.zeros(X.shape[1]);mean=np.array([np.median(raw[ii]) for ii in groupids])
 for iteration in range(12):
  def residual(p):
   sc=np.exp(X@p);I=raw*sc;s=expected_sigma(sc,mean,inv)
   return ((I-mean[inv])/np.sqrt(s*s+(.04*np.maximum(mean[inv],0))**2))[count0[inv]>1]
  fit=least_squares(residual,pscale,loss='soft_l1',f_scale=1.5,max_nfev=60)
  pscale=fit.x;sc=np.exp(X@pscale);I=raw*sc;s=expected_sigma(sc,mean,inv)
  mean=np.array([np.median(I[ii]) for ii in groupids])
 medsig=np.array([np.median(s[ii]) for ii in groupids])
 rejected=(count0[inv]>=3)&(mean[inv]>5*medsig[inv])&(abs(I-mean[inv])>np.maximum(8*s,.30*np.maximum(mean[inv],0)))
 # A bimodal even-sized group can put every observation far from the median.
 # In that case retain the group with inflated uncertainty, not an arbitrary mode.
 for ii in groupids:
  if rejected[ii].all():rejected[ii]=False
 w=(~rejected)/(s*s+(.04*np.maximum(mean[inv],0))**2)
 sw=np.bincount(inv,weights=w,minlength=ng)
 mean=np.bincount(inv,weights=w*I,minlength=ng)/sw
 res=(I-mean[inv])/np.sqrt(s*s+(.04*np.maximum(mean[inv],0))**2)
 count=np.bincount(inv,weights=(~rejected).astype(int),minlength=ng)
 scatter=np.bincount(inv,weights=w*(I-mean[inv])**2,minlength=ng)
 sem=np.sqrt(1/sw*np.maximum(1,scatter/np.maximum(count-1,1)))
 rint_before=np.sum(np.abs(I-mean[inv]))/np.sum(np.abs(I))
 rint=np.sum(np.abs(I-mean[inv])[~rejected])/np.sum(np.abs(I)[~rejected])
 strong=(mean[inv]>3*sem[inv])&~rejected
 rintstrong=np.sum(np.abs(I-mean[inv])[strong])/np.sum(np.abs(I)[strong])
 np.savetxt(OUT/'rejected_observations.csv',np.c_[data[rejected,:13],I[rejected],res[rejected]],delimiter=',',header='run,h,k,l,omega,x,y,module,Linv,P,edge,Iunscaled,sigmaunscaled,Iscaled,z',comments='',fmt='%.8g')
 flags=(a[:,10]!=0).astype(int)+2*(a[:,17]<2).astype(int)+4*(a[:,8]<=.025).astype(int)+8*shadow.astype(int)+16*(~np.isfinite(a[:,11])|~np.isfinite(a[:,12])|(a[:,12]<=0)).astype(int)
 flags[np.flatnonzero(valid)[rejected]]|=32
 np.savetxt(OUT/'observation_flags.csv',np.c_[np.arange(len(a)),a[:,:8],flags],delimiter=',',header='row_zero_based,run,h,k,l,omega,x,y,module,flags_1_edge_2_short_4_lowL_8_shadow_16_invalid_32_symmetry_outlier',comments='',fmt='%.8g')
 print('rejected observations',int(rejected.sum()),'Rint before rejection',rint_before,flush=True)
 cell=np.array(json.loads((OUT/'geometry.json').read_text())['cell']);d=1/np.linalg.norm(uniq/cell,axis=1)
 systematic=np.zeros(ng,bool)
 for axis in range(3):
  others=[j for j in range(3) if j!=axis]
  systematic|=(uniq[:,others]==0).all(axis=1)&(uniq[:,axis]%2==1)
 merged=np.c_[uniq,mean,sem,count,d,systematic.astype(int)]
 np.savetxt(OUT/'merged_all.csv',merged,delimiter=',',header='h,k,l,I,sigma,multiplicity,d,absent_P212121',comments='',fmt='%.8g')
 np.save(OUT/'merged.npy',merged)
 np.savetxt(OUT/'scaled_observations.csv',np.c_[data,I,s,res],delimiter=',',header='run,h,k,l,omega,x,y,module,Linv,P,edge,Iunscaled,sigmaunscaled,Iraw,varraw,counts,bgsum,nframes,wom,wpos,wom2,Iscaled,sigmascaled,normalized_residual',comments='',fmt='%.8g')
 with (OUT/'observed.hkl').open('w') as f:
  for row in merged[~systematic]:f.write('%4d%4d%4d%12.4f%12.4f\n'%tuple(row[:5]))
  f.write('   0   0   0      0.0000      0.0000\n')
 print('scale factors runs2..6,module1',np.exp(fit.x),flush=True)
 print('N raw/valid/unique',len(a),len(data),len(uniq),'Rint',rint,'Rint strong',rintstrong,flush=True)
 print('SYSTEMATIC ABSENCES',merged[systematic],flush=True)
 stats=dict(n_predicted=len(a),n_used=len(data),n_retained=len(data)-int(rejected.sum()),n_shadow_masked=int(shadow.sum()),n_geometry_filtered=int((~valid).sum()),counting_variance='Expected signal counts from equivalent-group median plus measured background variance; 4 percent group-intensity floor; group scatter inflation',n_unique=len(uniq),n_unique_allowed=int((~systematic).sum()),n_gt_2sigma=int(np.sum((mean>2*sem)&~systematic)),Rint=float(rint),Rint_strong=float(rintstrong),scale_parameters=fit.x.tolist(),scale_factors=np.exp(fit.x).tolist(),n_outlier_gt5sigma=int(np.sum(abs(res)>5)),outliers_omitted=True,n_rejected=int(rejected.sum()),Rint_before_rejection=float(rint_before),rejection_rule='at least 3 equivalents, median SNR>5, deviation from median greater than max(8sigma,30 percent median); symmetry-only, no structural model',median_I_sigma=float(np.median(mean[~systematic]/sem[~systematic])))
 # Geometrical completeness against all P212121-allowed mmm unique indices.
 ranges=[range(int(np.ceil(c/.70))+1) for c in cell]
 hall=np.array(list(itertools.product(*ranges)));norm=np.linalg.norm(hall/cell,axis=1)
 for cutoff in [.85,.80,.75,.72,.71,.70]:
  use=(norm>0)&(norm<=1/cutoff)
  for axis in range(3):
   others=[j for j in range(3) if j!=axis]
   use&=~((hall[:,others]==0).all(axis=1)&(hall[:,axis]%2==1))
  npossible=use.sum();nobs=((d>=cutoff)&~systematic).sum()
  stats[f'completeness_d{cutoff}']=dict(observed=int(nobs),possible=int(npossible),fraction=float(nobs/npossible))
 (OUT/'merging_statistics.json').write_text(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2),flush=True)
if __name__=='__main__':merge()
