exec(open('deg.py').read().split("Z=tangential_zeros")[0])
Z=tangential_zeros(dphi)
odd=all(any(np.linalg.norm(x+y)<1e-6 for y,_ in Z) for x,_ in Z)
print('zeros:',len(Z),' antipodal pairs:',odd)
eps=1e-2; D=dphi(X)
ev=np.unique(np.round([eps*l/Rp for _,l in Z],12)); edges=np.r_[-.05,ev,.05]
sgn=-1  # identity map has N=-1 in these Hopf coordinates; report oriented degree
for lo,hi in zip(edges,edges[1:]):
    t=(lo+hi)/2; R=(np.sqrt(3)/2)*np.cos(2*(np.pi/4+t))
    F=R*X+eps*D; print(f'  tau in ({lo:+.5f},{hi:+.5f}): N={sgn*degree(F,h):+.3f}  min|phi|={np.linalg.norm(F,axis=-1).min():.1e}')
# index of each zero event: sign det of spacetime Jacobian d(Phi)/d(tau,x) restricted, via local degree
