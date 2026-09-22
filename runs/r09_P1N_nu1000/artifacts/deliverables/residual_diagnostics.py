from crystal import *
from refine import metrics
import sys
name=sys.argv[1] if len(sys.argv)>1 else 'framework_aniso1'
d=np.load(OUT/'merged_equal.npz'); f=np.load(OUT/(name+'_fc.npz'));h=d['h'];I=d['I'];sig=d['sig'];q=d['q'];ic=f['ic']
rows=[]
for lo,hi in [(0,.1),(.1,.2),(.2,.3),(.3,.4),(.4,.5),(.5,.6),(.6,.7),(.7,.8),(.8,.9),(.9,1),(1,1.1),(1.1,1.2)]:
    m=(q>lo**2)&(q<=hi**2);pos=m&(I>2*sig);w=1/(sig**2+(.08*np.maximum(I,0))**2+.01)
    k=np.sum(w[m]*ic[m]*I[m])/np.sum(w[m]*ic[m]**2)
    row={'inv_d':[lo,hi],**metrics(I,sig,np.sqrt(ic),m),'relative_intensity_scale':float(k),'R1_rescaled':metrics(I,sig,np.sqrt(max(k,0)*ic),m)['R1_obs']}
    rows.append(row);print(row)
print('Strongest')
for j in np.argsort(I)[::-1][:40]: print(h[j],round(I[j],3),round(sig[j],3),round(ic[j],3),'d',round(1/np.sqrt(q[j]),3))
print('Largest standardized differences')
r=(I-ic)/np.sqrt(sig**2+(.08*np.maximum(I,0))**2+.01)
for j in np.argsort(abs(r))[::-1][:40]: print(h[j],round(I[j],3),round(sig[j],3),round(ic[j],3),'d',round(1/np.sqrt(q[j]),3),'z',round(r[j],2))
(OUT/(name+'_resolution_residuals.json')).write_text(json.dumps(rows,indent=2))
