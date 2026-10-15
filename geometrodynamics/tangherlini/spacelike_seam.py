"""Leading spacelike junction induced by a reciprocal-map scalar packet.

This is an asymptotic junction calculation, not a nonlinear spacetime or a
surface action. Unit-normal quantities are undefined on the null prefix.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import cumulative_trapezoid
from .characteristic_packet import geometry


def integrate(y,V):
    return cumulative_trapezoid(y,x=V,initial=0)


def deflected_seam(V,trace,b=1.,c=None,kappa5=1.):
    V,trace=np.asarray(V,dtype=float),np.asarray(trace,dtype=float)
    if b<=0 or kappa5<=0 or V.ndim!=1 or V.shape!=trace.shape or len(V)<3 or np.any(V<=0) or np.any(np.diff(V)<=0):
        raise ValueError('require positive b,kappa5 and increasing positive V with matching trace')
    c=b*b if c is None else c
    if c<=0: raise ValueError('require c>0')
    phip=np.gradient(trace/b**1.5,V,edge_order=2)
    S=kappa5*b/3*phip**2
    I0=integrate(S,V);B=integrate(V*S,V);C=integrate(V*V*S,V)
    # L=integral s log(V/s) S(s) ds; nonnegative in exact arithmetic.
    L=np.log(V/b)*B-integrate(V*np.log(V/b)*S,V)
    if np.min(L)<-1e-14*max(np.max(B),1e-300):
        raise ArithmeticError('negative positive-kernel integral')
    K=2/(b*np.e**2);Fh=4/np.e**2
    rin=-V*I0+B;rout=-B+C/V;difference=rin-rout
    x=(-V*I0+B+L)/(K*V)
    xp=-L/(K*V*V)
    xpp=(2*L-B)/(K*V**3)
    g=-c/V;gp=c/V**2;gpp=-2*c/V**3
    y=V/c*(difference/K-V*x)
    yp=xp/gp;ypp=xpp/gp-xp*gpp/gp**2
    H=-Fh*xp
    active=H>0
    # No regularization of the null segment: unit-frame formulas are NaN there.
    root=np.full_like(H,np.nan);root[active]=np.sqrt(H[active])
    length_minus=(1-B/(2*np.where(active,L,np.nan)))/(V*root)
    length_plus=(-1-B/(2*np.where(active,L,np.nan)))/(V*root)
    angle_minus=-(B+2*L)/(b*V*root)
    angle_plus=(-B+2*L)/(b*V*root)
    jump_length=-2/(V*root)
    jump_angle=2*V*root/b**2
    # Spacelike Israel convention: kappa5 S^a_b = [K^a_b]-delta^a_b[K].
    stress_length=-3*jump_angle/kappa5  # leading O(A), remainder O(A^3)
    stress_angle=-jump_length/kappa5  # only O(1/A); omitted remainder O(A)
    Qh=3*np.exp(-2)/b**2
    P=Qh/b**1.5*(V*trace-integrate(trace,V))
    Q=Qh/b**1.5*(g*trace-integrate(gp*trace,V))
    scalar_jump=2*xp/root*(P-Q/gp)  # O(A^2), after enforcing common O(A^3) trace
    # Codazzi at O(A): derivative of S^s_s plus 3 R'/b (S^s_s-S^angle_angle).
    # The product R' S^s_s is O(A^3), excluded consistently.
    Rprime=-B/V
    Hprime=2*b*(B-2*L)/V**3
    ds=-6/(kappa5*b*b)*(root+V*Hprime/(2*root))
    codazzi=ds-3*Rprime/b*stress_angle
    gauss=(length_plus+length_minus)/2*jump_angle+(angle_plus+angle_minus)/2*jump_length
    # Signed scalar delta-source for Box Phi=J delta(tau): J=-[n.Phi].
    return dict(V=V,trace=trace,S=S,B=B,L=L,x=x,xprime=xp,xsecond=xpp,
                y=y,yprime=yp,ysecond=ypp,g=g,gprime=gp,gsecond=gpp,
                H=H,rootH=root,active=active,incoming_delta_r=rin,outgoing_delta_r=rout,
                common_radius_shift=rin-K*V*x,common_radius_prime=Rprime,
                radius_match_residual=rin-K*V*x-(rout-K*g*y),
                temporal_match_residual=xp-gp*yp,
                K_length_minus=length_minus,K_length_plus=length_plus,
                K_angle_minus=angle_minus,K_angle_plus=angle_plus,
                jump_K_length=jump_length,jump_K_angle=jump_angle,
                S_length_leading=stress_length,S_angle_leading=stress_angle,
                densitized_S_angle=2/(kappa5*V),
                densitized_null_generator=-6*V/(kappa5*b**2),
                P_incoming= P,Q_outgoing=Q,scalar_normal_jump_leading=scalar_jump,
                scalar_delta_source_leading=-scalar_jump,
                codazzi_leading_residual=codazzi,gauss_leading_residual=gauss)


def field_derivative_source_control(packet,b=1.):
    """Independent one-sided derivatives of the two 2D fields at the seam.

    Only unboosted spherical packet data are accepted by this control.
    """
    U,V=packet['U'],packet['V']
    outU,outV=packet['out_U'],packet['out_V']
    r,_,_=geometry(U[:,None],V[None,:],b=b)
    phi=packet['field'][...,0]/r**1.5
    P=np.gradient(phi,U,axis=0,edge_order=2)[-1]
    r,_,_=geometry(outU[:,None],outV[None,:],b=b)
    phi=packet['outgoing'][...,0]/r**1.5
    Q=np.gradient(phi,outV,axis=1,edge_order=2)[:,0]
    seam=deflected_seam(V,packet['field'][-1,:,0],b=b)
    jump=2*seam['xprime']/seam['rootH']*(P-Q/seam['gprime'])
    return dict(P=P,Q=Q,scalar_normal_jump=jump)
