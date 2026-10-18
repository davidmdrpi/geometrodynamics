#!/usr/bin/env python
"""Geometric study, not a wave-driven collapse calculation. See docs/finite_neck_visualization.md."""
import argparse
import hashlib
import json
from pathlib import Path
import math
import numpy as np
from PIL import Image, ImageDraw, ImageFont

H, A, B0, BS, TC, END = .6, 1., .745071089, .045, 6.8, 10.
BG = (12, 22, 34)
FG = (222, 233, 241)
MUTED = (150, 174, 190)
TEAL = (85, 215, 204)
GOLD = (255, 195, 111)
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

def smooth(u):
    u = np.clip(u, 0., 1.)
    return u*u*(3-2*u)

def waist(t):
    return B0+(BS-B0)*smooth((t-2.6)/(TC-2.6))

def profile(x, t):
    """Fixed rims; cosh family blends into a rounded double cone before surgery."""
    b = waist(t)
    cosh = b*np.cosh(np.arccosh(A/b)*x/H)
    cone = np.sqrt(b*b+(A*A-b*b)*(x/H)**2)
    mix = smooth((.5-b)/.3)
    return (1-mix)*cosh+mix*cone

def deformation(x,t):
    """Prescribed signed display packet; neither transport PDE nor shoaling law."""
    c = -1.+.35*t
    return .08*np.exp(-((x-c)/.24)**2)*np.sin(22*(x-c))*np.maximum(0,1-(x/H)**2)

def branches(t,n=61):
    if t < TC:
        x = np.linspace(-H,H,n)
        r = profile(x,t)
        u = deformation(x,t)
        return [(x,r*(1+u),u,r)]
    # Imposed reconnection to two closed convex caps, then prescribed flattening.
    length = H*(.98-.96*smooth((t-TC)/(END-TC)))
    v = np.linspace(0,1,n)
    u = 1-v
    r = A*np.sqrt((u*u+.03*u)/1.03)
    return [(-H+length*v,r,np.zeros(n),r), (H-length*v,r,np.zeros(n),r)]

def phase(t):
    if t < 2.6: return '01  Incoming deformation', 'Fixed rims; prescribed signed packet'
    if t < 5.2: return '02  Narrowing waist', 'Contraction prescribed independently of the packet'
    if t < TC: return '03  Rounded double cone', 'Finite cutoff: b = 0.045 at the surgery frame'
    if t < 7.5: return '04  Imposed reconnection', 'Connected neck replaced by two closed caps'
    return '05  Cap retraction', 'Prescribed geometric relaxation; no energy-transfer claim'

def draw_frame(t, size=(960,720)):
    im=Image.new('RGB',size,BG);d=ImageDraw.Draw(im)
    font=lambda n,b=False: ImageFont.truetype(BOLD if b else FONT,n)
    d.text((28,20),'A finite neck: wave, narrowing, reconnection',font=font(25,True),fill=FG)
    d.text((28,58),'GEOMETRIC STUDY  /  packet and collapse prescribed; no coupled dynamics',font=font(14),fill=GOLD)
    title,subtitle=phase(t)
    d.text((28,100),title,font=font(20,True),fill=TEAL)
    d.text((28,131),subtitle,font=font(13),fill=MUTED)
    d.text((925,102),f't = {t:4.2f}',font=font(17),fill=FG,anchor='ra')
    d.text((30,170),'CUTAWAY SURFACE',font=font(12,True),fill=MUTED)
    def project(x,y,z):
        return (320+190*(1.35*x-.6*y),350+150*(-z+.25*x+.22*y))
    cells=[]
    angles=np.linspace(np.pi/3,5*np.pi/3,37)
    for x,r,u,base in branches(t):
        for i in range(len(x)-1):
            for j in range(len(angles)-1):
                vertices=[(x[k],r[k]*np.cos(angles[l]),r[k]*np.sin(angles[l]))
                          for k,l in [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]]
                depth=np.mean([.45*p[0]+p[1]+.2*p[2] for p in vertices])
                wave=(u[i]+u[i+1])/2/.08
                color=np.array(TEAL if wave>=0 else (217,125,179),float)
                strength=min(abs(wave),1)
                color=(1-strength)*np.array((98,136,157))+strength*color
                shade=.62+.28*(np.cos((angles[j]+angles[j+1])/2-.5)+1)/2
                color=tuple((color*shade).astype(int))
                cells.append((depth,[project(*p) for p in vertices],color))
    for _,poly,color in sorted(cells,key=lambda x:x[0]):d.polygon(poly,fill=color)
    # Exposed longitudinal edges make both the lumen and the closure visible.
    for x,r,u,base in branches(t):
        for phi in [angles[0],angles[-1]]:
            d.line([project(xx,rr*np.cos(phi),rr*np.sin(phi)) for xx,rr in zip(x,r)],fill=FG,width=2)
    for rim in [-H,H]:
        ring=[project(rim,np.cos(p),np.sin(p)) for p in np.linspace(0,2*np.pi,120)]
        d.line(ring,fill=GOLD,width=2)
    if t<TC:
        b=waist(t)
        ring=[project(0,b*np.cos(p),b*np.sin(p)) for p in angles]
        d.line(ring,fill=TEAL,width=3)
    d.text((32,517),'Front sector removed for visibility; rims stay fixed.',font=font(12),fill=MUTED)
    # Exact axial section and nominal profile, aligned to the finite rim separation.
    d.text((650,170),'AXIAL SECTION',font=font(12,True),fill=MUTED)
    proj=lambda x,r:(780+155*x/H,350-130*r)
    d.line([proj(-H,0),proj(H,0)],fill=(49,69,87),width=1)
    for rim in [-H,H]:d.line([proj(rim,-1.08),proj(rim,1.08)],fill=GOLD,width=2)
    for x,r,u,base in branches(t,121):
        for sign in [-1,1]:
            d.line([proj(xx,sign*rr) for xx,rr in zip(x,base)],fill=(74,93,110),width=1)
            d.line([proj(xx,sign*rr) for xx,rr in zip(x,r)],fill=TEAL,width=3)
    if t<TC:
        b=waist(t);d.line([proj(0,-b),proj(0,b)],fill=GOLD,width=2)
        d.text((790,350),f' b = {b:.3f}',font=font(13),fill=GOLD)
    else:d.text((710,350),'two closed caps',font=font(13),fill=GOLD)
    d.text((653,517),'Cosh family → rounded cones → caps',font=font(12),fill=MUTED)
    # Waist timeline: no fictitious continuous waist after the topology change.
    d.line([(28,552),(932,552)],fill=(47,68,86),width=1)
    d.text((30,570),'PRESCRIBED WAIST',font=font(12,True),fill=MUTED)
    xp=lambda q:225+q/END*690
    yp=lambda b:650-b/B0*60
    d.line([(225,650),(915,650)],fill=(74,93,110),width=1)
    d.line([(xp(q),yp(waist(q))) for q in np.linspace(0,TC,160)],fill=TEAL,width=3)
    d.line([(xp(TC),580),(xp(TC),657)],fill=GOLD,width=2)
    d.text((xp(TC)+8,580),'surgery',font=font(12),fill=GOLD)
    if t<TC:
        px,py=xp(t),yp(waist(t));d.ellipse((px-5,py-5,px+5,py+5),fill=FG)
    d.line([(xp(t),625),(xp(t),657)],fill=FG,width=1)
    d.text((30,601),'Signed wave: inward / outward',font=font(12),fill=MUTED)
    d.text((30,623),'Display deformation ≤ 8%',font=font(12),fill=MUTED)
    d.text((28,684),'No inevitable pinch-off assumed. No singular radius is evaluated. No released energy is assigned.',font=font(13),fill=MUTED)
    return im

def verify():
    for t in np.linspace(0,TC,100,endpoint=False):
        assert np.max(abs(profile(np.array([-H,H]),t)-A))<1e-12
        assert abs(profile(0,t)-waist(t))<1e-12
        for x,r,u,base in branches(t):assert np.all(np.isfinite(r)) and np.min(r)>0
    for t in np.linspace(TC,END,20):
        parts = branches(t)
        assert parts[0][0][-1] < parts[1][0][-1]
        for x,r,u,base in parts:
            assert r[-1]==0 and r[0]==A and np.all(np.isfinite(r))
    print('Geometry checks: fixed rims, nominal waist, positive pre-surgery radius, closed post-surgery caps')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=Path('docs/images/finite_neck'))
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True);verify()
    times=np.unique(np.r_[np.linspace(0,END,51),TC-.03,TC])
    frames=[draw_frame(float(t)).resize((720,540)) for t in times]
    palette=Image.new('RGB',(720,540*6))
    for i,j in enumerate([0,13,25,32,36,50]):palette.paste(frames[j],(0,540*i))
    palette=palette.quantize(colors=64)
    quant=[f.quantize(palette=palette,dither=Image.Dither.NONE) for f in frames]
    duration=[100]*len(times)
    duration[np.flatnonzero(times==TC)[0]]=1100;duration[-1]=1800
    quant[0].save(args.output/'finite_neck.gif',save_all=True,append_images=quant[1:],duration=duration,loop=0,optimize=True)
    panels=[draw_frame(t) for t in [1.7,6.4,7.1]]
    sheet=Image.new('RGB',(1440,360),BG)
    for i,p in enumerate(panels):sheet.paste(p.resize((480,360)),(480*i,0))
    sheet.save(args.output/'sequence.png',optimize=True)
    draw_frame(6.4).save(args.output/'closeup.png',optimize=True)
    files=[args.output/x for x in ['finite_neck.gif','sequence.png','closeup.png']]
    (args.output/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2,sort_keys=True)+'\n')
    print({p.name:p.stat().st_size for p in files})

if __name__=='__main__':main()
