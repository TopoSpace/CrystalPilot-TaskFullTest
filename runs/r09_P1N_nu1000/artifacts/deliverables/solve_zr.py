from crystal import *
from scipy.optimize import least_squares
import os

def coords(p,typ):
    u,v,w=p[:3]
    if typ==0:return [[.5+u,0,w],[.5+v,2*v,0]]
    if typ==1:return [[.5+v,2*v,w],[.5+u,0,0]]
    return [[.5+u,0,0],[.5+v,2*v,0],[.5,0,w]]

def op_reps(x):
    eq=np.einsum('sij,j->si',R,np.array(x))%1
    _,idx=np.unique(np.round(eq,7),axis=0,return_index=True)
    return R[idx]

def run():
    d=np.load(OUT/'merged_equal.npz'); hall=d['h'];Iall=d['I'];sigall=d['sig'];qall=d['q']
    mask=(qall>.16**2)&(qall<.72**2)&(Iall>2*sigall)
    h=hall[mask]; I=Iall[mask];sig=sigall[mask];q=qall[mask]
    fo=np.sqrt(I);sf=sig/(2*fo); weight=1/np.sqrt(sf**2+(.12*fo)**2+.01)
    ff=formfactor('Zr',q/4)
    best=None
    for typ in [0,1,2]:
        reps=[op_reps(x) for x in coords([.045,.037,.11],typ)]
        print('Type',typ,'multiplicities',[len(x) for x in reps],flush=True)
        hr=[np.einsum('ni,sij->nsj',h,r) for r in reps]
        def fc(p,allh=False):
            hh=hr if not allh else [np.einsum('ni,sij->nsj',hall,r) for r in reps]
            qq=q if not allh else qall
            f=ff if not allh else formfactor('Zr',qall/4)
            sums=0
            for x,hrx in zip(coords(p,typ),hh):sums+=np.cos(2*np.pi*np.einsum('nsj,j->ns',hrx,x)).sum(axis=1)
            return np.exp(p[4])*f*np.exp(-p[3]*qq/4)*sums
        def fun(p):return (abs(fc(p))-fo)*weight
        rng=np.random.default_rng(374+typ)
        for j in range(14):
            if j==0:p0=[.045,.037,.11,3.,-3.]
            else:p0=[rng.uniform(.025,.075),rng.uniform(.02,.055),rng.uniform(.085,.165),3.,-3.]
            opt=least_squares(fun,p0,bounds=([.012,.008,.045,.1,-8],[.1,.075,.23,16,0]),max_nfev=180,ftol=1e-7)
            score=np.mean(fun(opt.x)**2)
            if best is None or score<best['score']:
                ffull=fc(opt.x,True)
                best={'score':score,'type':typ,'p':opt.x.tolist(),'sites':coords(opt.x,typ),'B':float(opt.x[3]),'amplitude_scale':float(np.exp(opt.x[4]))}
                (OUT/'zr_solution.json').write_text(json.dumps(best,indent=2))
                np.savez_compressed(OUT/'zr_phases.npz',h=hall,fc=ffull)
                print('NEW BEST',best,flush=True)
            print('trial',typ,j,'score',round(score,4),'params',np.round(opt.x,6),flush=True)
    print('FINAL',best,flush=True)
    # Difference map on physically useful resolution only.
    p=np.load(OUT/'zr_phases.npz');fc=p['fc'];scale=best['amplitude_scale']
    m=qall<1.05**2
    fobs=np.sqrt(np.maximum(Iall,0))
    for title,coeff in [('zr_fo',(fobs*np.sign(fc))/scale),('zr_diff',(fobs*np.sign(fc)-fc)/scale)]:
        rho=map_fft(hall[m],coeff[m],shape=(240,240,120))
        save_peaks(peaks(rho,n=150,min_dist=.65),title+'_peaks.csv')
        np.save(OUT/(title+'_map.npy'),rho.astype(np.float32))
        print(title,'peaks')
        for x,v in peaks(rho,n=35,min_dist=.7):print(np.round(x,6),round(v,3))
if __name__=='__main__':run()
