"""Predict reflections from the Ewald equation and integrate raw photon counts.
No crystallographic library is used. All negative net intensities are retained.
"""
import json,time,itertools,os
import numpy as np
from raw_frames import OUT,INP,read_frame
from geometry import LAMBDA,ub_runs,detector_basis
from geometry_probe import rotate_z

RADIUS=float(os.environ.get('ALANINE_RADIUS','4.5'))
ROCK_MIN=float(os.environ.get('ALANINE_ROCK_MIN','1.25'))
ROCK_HALF=10.0
DMIN=.70

def rock_width(gnorm,lorentz):
 return np.maximum(ROCK_MIN,np.minimum(10,.25+3.5*np.sqrt(.07**2+(.10*gnorm/np.maximum(lorentz,.01))**2)))

def predictions(run,par,meta):
 ub=ub_runs(par)[run-1]
 cell=par[:3];limits=np.ceil(cell/DMIN).astype(int)
 h=np.array(list(itertools.product(*[range(-n,n+1) for n in limits])),int)
 v=h@ub.T;norm2=np.sum(v*v,axis=1)
 good=(norm2>0)&(norm2<(LAMBDA/DMIN)**2)
 h=h[good];v=v[good];norm2=norm2[good]
 rho=np.hypot(v[:,0],v[:,1]);cosang=-norm2/2/np.maximum(rho,1e-12)
 good=np.abs(cosang)<=1;h=h[good];v=v[good];cosang=cosang[good]
 phase=np.arctan2(v[:,1],v[:,0]);delta=np.arccos(cosang)
 oms=np.rad2deg(np.r_[phase+delta,phase-delta]);h=np.r_[h,h];v=np.r_[v,v]
 oms=(oms+180)%360-180
 th=meta[0]['theta'];omin=min(m['omega0'] for m in meta);omax=max(m['omega1'] for m in meta)
 good=(oms>omin-ROCK_HALF)&(oms<omax+ROCK_HALF)
 h=h[good];v=v[good];oms=oms[good]
 g=rotate_z(v,-oms);s=g+[1,0,0]
 sd=rotate_z(s,np.full(len(s),th))
 mod,normal,ex,ey=detector_basis(par)
 rows=[]
 for j in range(2):
  t=mod[j,0]/(sd@normal[j]);dx=t*(sd@ex[j])/.1;dy=t*(sd@ey[j])/.1
  x=mod[j,1]+dx;y=mod[j,2]+dy
  good=(x>RADIUS+1+j*415)&(x<(384 if j==0 else 799)-RADIUS-1)&(y>RADIUS+1)&(y<774-RADIUS-1)&(t>0)
  for i in np.where(good)[0]:
   lorentz=abs(g[i,1]);pol=(1+s[i,0]**2)/2
   half=rock_width(np.linalg.norm(g[i]),lorentz)
   rows.append([run,*h[i],oms[i],x[i],y[i],j,lorentz,pol,1 if oms[i]-half<omin or oms[i]+half>omax else 0])
 return np.array(rows)

def aperture(im,x,y,radius=RADIUS):
 ix=int(round(x));iy=int(round(y));inner=max(radius+2,7);outer=max(radius+5,10);r=int(np.ceil(outer))
 x0=max(0,ix-r);x1=min(800,ix+r+1);y0=max(0,iy-r);y1=min(775,iy+r+1)
 box=im[y0:y1,x0:x1];yy,xx=np.mgrid[y0:y1,x0:x1]
 rr=(xx-x)**2+(yy-y)**2
 sig=(rr<=radius**2)&(box>=0)
 bg=(rr>=inner**2)&(rr<=outer**2)&(box>=0)&((xx<383)|(xx>416))
 vals=box[bg]
 # Sparse Poisson backgrounds: use a mean, never a zero-biased median.
 for _ in range(2):
  if not len(vals):return None
  mu=vals.mean();vals=vals[vals<=max(5,mu+5*np.sqrt(mu+1))]
 n=sig.sum();raw=box[sig].sum();back=vals.mean();net=raw-n*back
 var=max(raw,0)+n*n*max(back,.01)/len(vals)
 return net,var,raw,back,n

def integrate_run(run):
 par=np.array(json.loads((OUT/'geometry.json').read_text())['parameters'])
 meta=json.loads((OUT/'frame_metadata.json').read_text());meta=[m for m in meta if int(m['file'].split('_')[-2])==run]
 pred=predictions(run,par,meta);acc=np.zeros((len(pred),8));profiles=[]
 halfwidth=rock_width(np.linalg.norm(pred[:,1:4]/par[:3],axis=1)*LAMBDA,pred[:,8])
 for m in meta:
  om=(m['omega0']+m['omega1'])/2
  use=np.where(np.abs(pred[:,4]-om)<=halfwidth+.001)[0]
  if not len(use):continue
  im,_=read_frame(INP/'frames'/m['file'])
  for i in use:
   z=aperture(im,pred[i,5],pred[i,6])
   if z is None:continue
   net,var,raw,bg,npix=z
   acc[i,0]+=net;acc[i,1]+=var;acc[i,2]+=raw;acc[i,3]+=bg;acc[i,4]+=1
   # intensity-weighted rocking centroid diagnostic only (not used in net intensity)
   acc[i,5]+=max(net,0)*om;acc[i,6]+=max(net,0);acc[i,7]+=max(net,0)*om*om
   profiles.append([i,int(m['file'].split('_')[-1].split('.')[0]),om,net,var,bg])
 corr=pred[:,8]/pred[:,9]
 data=np.c_[pred,acc[:,0]*corr,np.sqrt(acc[:,1])*corr,acc]
 np.save(OUT/f'integrated_run{run}.npy',data)
 np.save(OUT/f'profiles_run{run}.npy',np.array(profiles))
 print('integrated run',run,'predictions',len(pred),'positive>3sigma',int(np.sum(data[:,11]>3*data[:,12])),flush=True)
 return len(pred)

def main():
 from concurrent.futures import ProcessPoolExecutor
 t=time.time()
 with ProcessPoolExecutor(max_workers=4) as pool:ns=list(pool.map(integrate_run,range(1,7)))
 data=np.concatenate([np.load(OUT/f'integrated_run{r}.npy') for r in range(1,7)])
 header='run,h,k,l,omega,x,y,module,Linv,polarization,scan_edge,I,sigma,Iraw_net,var_raw,raw_counts,bg_sum,nframes,weighted_omega,sum_pos,weighted_omega2'
 np.savetxt(OUT/'unmerged_integrated.csv',data,delimiter=',',header=header,comments='',fmt='%.8g')
 print('integrated total',len(data),'seconds',round(time.time()-t,2),flush=True)
 (OUT/'integration_settings.json').write_text(json.dumps(dict(aperture_radius_pixels=RADIUS,rock_halfwidth_degrees=f'max({ROCK_MIN},min(10,0.25+3.5*sqrt(0.07**2+(0.10*norm(g)/abs(g_y))**2)))',d_min_angstrom=DMIN,background_inner_radius=max(RADIUS+2,7),background_outer_radius=max(RADIUS+5,10),lorentz='multiply abs(g_y) for omega rotation about z',polarization='divide (1+s_x**2)/2',intensity_threshold='none in prediction or integration'),indent=2))
if __name__=='__main__':main()
