"""Constrained spherical perturbations of the #331 canonical scalar/EGB shell.

K and M below are an auxiliary self-adjoint spectral pencil. Its M is NOT the
physical kinetic metric. The physical metric is 6*C*R**4*(K+e*M), fixed by the
bulk plus Myers/surface symplectic form; see docs/egb_coupled_radial.md.
"""
from dataclasses import dataclass, field
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp, quad
from scipy.linalg import eigvals, eigvalsh
from scipy.optimize import brentq
from .egb_shell_action import ScalarShellAction


@dataclass(frozen=True)
class CoupledRadial:
    action: ScalarShellAction = field(default_factory=ScalarShellAction)

    def __post_init__(self):
        if self.action.bulk_lambda != 0:
            raise ValueError('this reduction is for the canonical bulk lambda=0')
        if self.D >= 0:
            raise ValueError('this implementation uses the nondegenerate D<0 branch')

    @property
    def shell(self): return self.action.shell
    @property
    def R(self): return self.shell.radius
    @property
    def a(self): return self.shell.mouth_angle
    @property
    def r(self): return self.shell.r0
    @property
    def C(self): return self.shell.bulk_C
    @property
    def D(self): return float(self.shell.linear_data()['D'])
    @property
    def W(self): return sum(self.shell.equilibrium)
    @property
    def x(self): return math.sin(self.a)**2
    @property
    def e(self): return 2/self.C
    @property
    def d(self): return 4/(math.cos(self.a)*self.D)
    @property
    def b(self): return self.R*self.W/(6*self.C*self.shell.eta)
    @property
    def j(self): return self.R*self.D/(2*self.C*self.r)
    @property
    def mid(self): return math.pi/2

    def shoot(self, eigenvalue, parity='even', tolerance=2e-11):
        """Eigenvalue is omega**2*R**2, not spatial Laplacian eigenvalue."""
        if parity not in ('even','odd'): raise ValueError('parity must be even or odd')
        initial=[1.,0.] if parity=='even' else [0.,1.]
        result=solve_ivp(lambda u,y:[y[1],-3/math.tan(u)*y[1]-(eigenvalue+4+self.e)*y[0]],
                         (self.mid,self.a),initial,method='DOP853',
                         rtol=tolerance,atol=tolerance/20,dense_output=True)
        if not result.success: raise RuntimeError(result.message)
        return result

    def determinant(self, eigenvalue, parity='even', tolerance=2e-11):
        f,fp=self.shoot(eigenvalue,parity,tolerance).y[:,-1]
        if parity=='odd': return float(f)
        h=eigenvalue+self.e; g=self.x*eigenvalue+self.d
        # Pole-free coupled determinant; never divide by g during root finding.
        return float(g*(fp+math.tan(self.a)*f+self.b*h*f)-self.j*h*f)

    def root(self, bracket, parity='even', tolerance=2e-11):
        return brentq(lambda z:self.determinant(z,parity,tolerance),*bracket,xtol=2e-12)

    def mode(self, eigenvalue, parity='even', tolerance=2e-11):
        solution=self.shoot(eigenvalue,parity,tolerance)
        if parity=='even':
            # Unit fractional areal shell displacement z, hence physical norms
            # are comparable without an arbitrary eigenvector normalization.
            h=eigenvalue+self.e; g=self.x*eigenvalue+self.d
            if abs(h)<1e-9: raise ValueError('exceptional h=0 requires separate gauge analysis')
            scale=(g/h)/solution.y[0,-1]
            z=1.
        else:
            # Odd modes do not move the shell; unit half-bulk weighted L2 norm.
            scale=1/math.sqrt(quad(lambda u:math.sin(u)**3*solution.sol(u)[0]**2,self.a,self.mid)[0])
            z=0.
        return lambda u:scale*solution.sol(u),z

    def physical_norm(self, eigenvalue, parity='even'):
        profile,z=self.mode(eigenvalue,parity);fa=float(profile(self.a)[0])
        integral=quad(lambda u:math.sin(u)**3*(profile(u)[1]**2-4*profile(u)[0]**2),
                      self.a,self.mid,epsabs=1e-10,epsrel=2e-10)[0]
        return float(6*self.C*self.R**4*(integral-math.sin(self.a)**3*math.tan(self.a)*fa**2)
                     +3*self.D*self.r**4*z*z)

    def spectral_norm(self, eigenvalue):
        profile,_=self.mode(eigenvalue);fa=float(profile(self.a)[0]);h=eigenvalue+self.e
        g=self.x*eigenvalue+self.d
        Fprime=-self.b+self.j*(self.d-self.x*self.e)/g**2
        I=quad(lambda u:math.sin(u)**3*profile(u)[0]**2,self.a,self.mid,epsabs=1e-11)[0]
        return float(6*self.C*self.R**4*h*(I-math.sin(self.a)**3*Fprime*fa*fa))

    def symplectic_norm(self,eigenvalue,gauge_power=1):
        """Direct ADM+scalar+Myers norm in a fixed-shell, zero-shift gauge.

        The arbitrary radial gauge profile E vanishes at the equator. Different
        E must give the same total, although individual bulk/surface parts vary.
        """
        if gauge_power not in (1,3,5): raise ValueError('use a smooth odd reflection gauge')
        profile,z=self.mode(eigenvalue);fa=float(profile(self.a)[0]);h=eigenvalue+self.e
        xi=(z-fa)*math.tan(self.a);width=self.mid-self.a;p=gauge_power
        E=lambda u:((self.mid-u)/width)**p
        Ep=lambda u:-p*(self.mid-u)**(p-1)/width**p
        J=lambda u:-(self.mid-u)**(p+1)/((p+1)*width**p)
        def integrand(u):
            f=profile(u)[0];ell=f+xi*Ep(u);psi=f+xi*E(u)/math.tan(u)
            scalar=3*self.C*((eigenvalue+4+self.e)*f+xi*(Ep(u)+3*E(u)/math.tan(u)))*(f-xi*J(u))
            gravity=-6*self.C*(ell*psi+psi*psi)
            return self.R**4*math.sin(u)**3*(scalar+gravity)
        bulk=2*quad(integrand,self.a,self.mid,epsabs=1e-9,epsrel=2e-10)[0]
        surface=self.r**3*self.R**2*self.W*(3*z+h*fa/self.shell.eta)*(fa-xi*J(self.a))
        myers=48*self.shell.alpha*self.r**3/(self.R*math.tan(self.a))*z*z
        return dict(bulk=float(bulk),surface_scalar=float(surface),myers=float(myers),total=float(bulk+surface+myers))

    def charge_residual(self,eigenvalue):
        profile,z=self.mode(eigenvalue);fa=float(profile(self.a)[0]);h=eigenvalue+self.e
        q=math.sqrt(self.action.q2);wa=math.sin(self.a)**3;xi=(z-fa)*math.tan(self.a)
        I=quad(lambda u:math.sin(u)**3*profile(u)[0],self.a,self.mid,epsabs=1e-11)[0]
        bulk=2*q*self.R**4*((eigenvalue+4+self.e)*I-wa*xi)
        surface=self.r**3/q*self.W*(3*z+h*fa/self.shell.eta)
        return float(bulk+surface)

    def junction_residuals(self,eigenvalue):
        profile,z=self.mode(eigenvalue);fa,fp=profile(self.a);h=eigenvalue+self.e
        xi=(z-fa)*math.tan(self.a);Q0=1/(self.R*math.tan(self.a))
        dQ=(fp-fa/math.tan(self.a)-xi/math.sin(self.a)**2)/self.R
        dk0=-eigenvalue*xi/self.R-2*fp/(self.C*self.R)
        ds=-6*self.C*dQ+48*self.shell.alpha*Q0*z/self.r**2
        dp=4*dQ+2*self.C*dk0-16*self.shell.alpha*Q0*eigenvalue*z/self.R**2
        material=self.W*h*fa
        current=material/self.shell.eta+3*self.W*z-6*self.C/self.R*(-fp+xi)
        return dict(density=float(ds-material/self.shell.eta),pressure=float(dp-material),scalar_current_integrated=float(current))

    def pencil(self,elements):
        """Independent linear finite-element discretization, no shooting input."""
        if elements<4: raise ValueError('need at least four elements')
        nodes=np.linspace(self.a,self.mid,elements+1)
        M=np.zeros((elements+2,elements+2));K=M.copy();qs,ws=leggauss(5)
        for i,(u0,u1) in enumerate(zip(nodes[:-1],nodes[1:])):
            dx=u1-u0;u=(u0+u1)/2+dx/2*qs;w=ws*dx/2*np.sin(u)**3
            shape=np.array([(1-qs)/2,(1+qs)/2]);ds=np.array([-1.,1.])/dx
            mass=np.einsum('iq,jq,q->ij',shape,shape,w)
            K[i:i+2,i:i+2]+=np.outer(ds,ds)*sum(w)-(4+self.e)*mass
            M[i:i+2,i:i+2]+=mass
        wa=math.sin(self.a)**3;My=-wa*self.j/(self.d-self.x*self.e)
        M[0,0]+=wa*self.b;M[-1,-1]=My
        K[0,0]+=wa*(-math.tan(self.a)-self.b*self.e+self.j/self.x)
        K[0,-1]=K[-1,0]=wa*self.j/self.x;K[-1,-1]=My*(-self.d/self.x)
        return K,M

    def finite_element_spectrum(self,elements):
        K,M=self.pencil(elements);values=eigvals(K,M);kinetic=eigvalsh(K+self.e*M)
        return dict(elements=elements,low_eigenvalues=sorted(float(v.real) for v in values if abs(v.imag)<1e-7)[:8],
                    maximum_imaginary_part=float(max(abs(values.imag))),
                    physical_negative_directions=int(sum(kinetic < -1e-9)),
                    auxiliary_negative_directions=int(sum(eigvalsh(M)<-1e-9)))


def egb_tensor(sectionals,alpha,mixed=None):
    """Independent orthonormal Riemann contraction, signature -++++.

    sectionals maps (a,b) to R_abab. No reduced field equations are used.
    """
    signature=np.array([-1.,1.,1.,1.,1.]);g=np.diag(signature);curv=np.zeros((5,5,5,5))
    for (a,b),v in sectionals.items():
        curv[a,b,a,b]=curv[b,a,b,a]=v
        curv[a,b,b,a]=curv[b,a,a,b]=-v
    for (a,b,c,d),value in (mixed or {}).items():
        for i,j,k,l,sign in [(a,b,c,d,1),(b,a,c,d,-1),(a,b,d,c,-1),(b,a,d,c,1),
                             (c,d,a,b,1),(d,c,a,b,-1),(c,d,b,a,-1),(d,c,b,a,1)]:
            curv[i,j,k,l]=sign*value
    ric=np.einsum('a,abad->bd',signature,curv);scalar=float(np.sum(signature*np.diag(ric)))
    riem2=np.einsum('a,b,c,d,abcd,abcd',signature,signature,signature,signature,curv,curv)
    ric2=np.einsum('a,b,ab,ab',signature,signature,ric,ric);gb=riem2-4*ric2+scalar*scalar
    H=2*(scalar*ric-2*np.einsum('c,ac,bc->ab',signature,ric,ric)
         -2*np.einsum('c,d,cd,acbd->ab',signature,signature,ric,curv)
         +np.einsum('c,d,e,acde,bcde->ab',signature,signature,signature,curv,curv))-.5*g*gb
    return ric-.5*g*scalar+alpha*H
