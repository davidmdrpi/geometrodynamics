#!/usr/bin/env python
"""Wave-contact visualization: prescribed packet, contact-triggered regularized reconnection."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from neck_wave_contact import ContactNeck, Parameters, H, END, X, smooth

BG=(12,22,34);FG=(222,233,241);MUTED=(150,174,190)
TEAL=(85,215,204);GOLD=(255,195,111);PINK=(227,133,182)
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

def project(x,y,z):
    return (310+195*(1.35*x-.6*y),330+150*(-z+.25*x+.22*y))

def arrow(d,points,color=GOLD,width=2):
    d.line(points,fill=color,width=width)
    v=np.array(points[-1])-np.array(points[-2]);v/=max(np.linalg.norm(v),1e-12)
    p=np.array(points[-1]);n=np.array([-v[1],v[0]])
    d.polygon([tuple(p),tuple(p-7*v+4*n),tuple(p-7*v-4*n)],fill=color)

def surface(parts):
    im=Image.new('RGB',(960,720),BG);d=ImageDraw.Draw(im);cells=[]
    phis=np.linspace(np.pi/3,5*np.pi/3,33)
    for x,r,u,_ in parts:
        for i in range(len(x)-1):
            for j in range(len(phis)-1):
                v=[(x[k],r[k]*np.cos(phis[l]),r[k]*np.sin(phis[l])) for k,l in [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]]
                dep=sum(.45*p[0]+p[1]+.2*p[2] for p in v)/4
                wave=(u[i]+u[i+1])/.36;f=min(abs(wave),1)
                col=np.array(TEAL if wave>=0 else PINK)*f+np.array((83,126,151))*(1-f)
                shade=.64+.26*(1+np.cos((phis[j]+phis[j+1])/2-.5))/2
                cells.append((dep,[project(*p) for p in v],tuple((col*shade).astype(int))))
    for _,p,c in sorted(cells,key=lambda v:v[0]):d.polygon(p,fill=c)
    for x,r,_,_ in parts:
        for phi in [phis[0],phis[-1]]:
            d.line([project(xx,rr*np.cos(phi),rr*np.sin(phi)) for xx,rr in zip(x,r)],fill=FG,width=2)
    return im

def draw_frame(m,t,curves):
    parts=m.branches(t,81);im=surface(parts)
    dt=-1 if m.tc is None else t-m.tc
    if 0<=dt<.5:
        # Explicit finite-duration rendering blend between two topologies.
        im=Image.blend(surface(m.branches(m.tc-1e-8,81)),im,float(smooth(dt/.5)))
    d=ImageDraw.Draw(im);font=lambda n,b=False:ImageFont.truetype(BOLD if b else FONT,n)
    d.text((25,18),'The deformation meets itself in a narrow neck',font=font(25,True),fill=FG)
    d.text((25,55),'FIXED RESTING WAIST  /  contact comes from the prescribed wave',font=font(14),fill=GOLD)
    title,note=m.phase(t);d.text((25,96),title,font=font(21,True),fill=TEAL)
    d.text((25,128),note,font=font(13),fill=MUTED)
    d.text((930,99),f't = {t:.2f}',font=font(17),fill=FG,anchor='ra')
    d.text((28,166),'CUTAWAY: signed deformation on a fixed narrow bridge',font=font(12),fill=MUTED)
    for rim in [-H,H]:d.line([project(rim,np.cos(p),np.sin(p)) for p in np.linspace(0,2*np.pi,100)],fill=GOLD,width=2)
    if m.tc is not None and -.5<dt<.7:
        rad=m.par.thickness+.008
        d.line([project(m.xc,rad*np.cos(p),rad*np.sin(p)) for p in np.linspace(0,2*np.pi,100)],fill=GOLD,width=5)
    if dt>=.5:
        # Arrows trace the imposed cap profiles. They are not velocity measurements.
        for x,r,_,_ in parts:
            j=int(15+35*((dt*.6)%1))
            pts=[project(x[k],r[k]*np.cos(1.3),r[k]*np.sin(1.3)) for k in [j+4,j+2,j]]
            arrow(d,pts)
    d.text((28,503),'No independent shrinking-waist schedule.',font=font(13),fill=FG)
    # Magnified meridional section, including finite-thickness facing surfaces.
    d.text((645,166),'CONTACT CLOSE-UP',font=font(13,True),fill=MUTED)
    center=m.xc if m.xc is not None else 0.
    halfwidth=.23;vertical=.26
    if dt>=0:
        tips=[part[0][-1] for part in parts]
        center=sum(tips)/2
        halfwidth=max(.23,(tips[1]-tips[0])/2+.13)
        vertical=.26+.7*m.progress(t)
    zlo,zhi=center-halfwidth,center+halfwidth
    box=(635,188,935,492);zoom=Image.new('RGB',(300,304),(16,30,43));zd=ImageDraw.Draw(zoom)
    proj=lambda x,y:((x-zlo)/(zhi-zlo)*300,152-y/vertical*145)
    z=np.linspace(zlo,zhi,501);eps=m.par.thickness
    if m.tc is None or t<m.tc:
        r=m.radius(z,t)
        for s in [-1,1]:
            poly=[proj(xx,s*(rr+eps)) for xx,rr in zip(z,r)]+[proj(xx,s*max(rr-eps,0)) for xx,rr in zip(z[::-1],r[::-1])]
            zd.polygon(poly,fill=(111,78,115));zd.line([proj(xx,s*rr) for xx,rr in zip(z,r)],fill=PINK,width=2)
            zd.line([proj(xx,s*rr) for xx,rr in zip(z,m.rest(z))],fill=(59,85,104),width=1)
    else:
        for x,r,_,_ in parts:
            for s in [-1,1]:
                pts=[proj(xx,s*rr) for xx,rr in zip(x,r)]
                zd.line(pts,fill=TEAL,width=5)
                # Curved motion guides near the rounded shoulder.
                if dt<1.5:
                    j=len(x)-12;px,py=proj(x[j],s*r[j])
                    if 10<px<290 and 15<py<289:
                        zd.arc((px-11,py-11,px+11,py+11),40,290,fill=GOLD,width=2)
        if dt<.5:
            # Contact is a finite diffuse region, not a rendered point crossing.
            a=float(1-smooth(dt/.5));layer=Image.new('RGB',(300,304),(16,30,43));ld=ImageDraw.Draw(layer)
            ld.ellipse((118,119,182,185),fill=GOLD)
            zoom=Image.blend(zoom,layer,.5*a);zd=ImageDraw.Draw(zoom)
    im.paste(zoom,(635,188));d=ImageDraw.Draw(im)
    d.rectangle(box,outline=(54,80,99),width=1)
    d.text((645,502),'Contact band → rounded rims (zoom follows)',font=font(12),fill=MUTED)
    # Contact-gap history includes no-wave and weak-wave controls.
    d.line([(25,545),(935,545)],fill=(47,68,86),width=1)
    d.text((25,560),'FACING-SURFACE GAP',font=font(12,True),fill=MUTED)
    px=lambda q:260+q/END*660;py=lambda gap:650-gap/.28*70
    d.line([(260,650),(920,650)],fill=GOLD,width=1)
    colors=[TEAL,(153,177,198),(77,98,118)]
    for curve,col in zip(curves,colors):d.line([(px(q),py(g)) for q,g in curve],fill=col,width=2)
    d.text((25,588),'Teal: wave reaches contact',font=font(12),fill=TEAL)
    d.text((25,610),'Gray: weak wave / wave off',font=font(12),fill=MUTED)
    if m.tc is not None:
        d.line([(px(m.tc),570),(px(m.tc),655)],fill=GOLD,width=1)
        d.text((px(m.tc)+8,561),f'contact t={m.tc:.3f}',font=font(12),fill=GOLD)
    d.line([(px(t),630),(px(t),655)],fill=FG,width=2)
    d.text((25,683),'Kinematic contact rule and illustrative rim motion. No fluid solver, energy release, or smooth-topology theorem.',font=font(12),fill=MUTED)
    return im

def verify():
    for amp in [0,.08]:assert ContactNeck(Parameters(amplitude=amp)).tc is None
    m=ContactNeck();assert m.tc is not None and abs(m.gap(m.tc))<1e-10
    assert m.gap(m.tc-.001)>0
    for t in np.linspace(0,m.tc-1e-8,50):
        assert min(m.radius(X,t))>0
        assert max(abs(m.radius(np.array([-H,H]),t)-1))<1e-12
    for t in np.linspace(m.tc,END,20):
        a,b=m.branches(t);assert a[0][-1]<b[0][-1]
        for x,r,_,_ in [a,b]:assert r[0]==1 and r[-1]==0 and np.all(r>=0)
    print('Verified: wave-off/weak-wave nulls, contact bracket, fixed rims, positive pre-contact neck and separated rounded caps')
    return m

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=Path('docs/images/finite_neck'));args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    m=verify();models=[m,ContactNeck(Parameters(amplitude=.08)),ContactNeck(Parameters(amplitude=0))]
    curves=[[(t,max(mm.gap(t),0)) for t in np.linspace(0,mm.tc if mm.tc is not None else END,161)] for mm in models]
    times=np.unique(np.r_[np.linspace(0,END,53),m.tc-.03,m.tc,m.tc+.25,m.tc+.5])
    frames=[draw_frame(m,float(t),curves).resize((720,540)) for t in times]
    palette=Image.new('RGB',(720,540*6))
    for i,j in enumerate(np.linspace(0,len(times)-1,6).astype(int)):palette.paste(frames[j],(0,540*i))
    palette=palette.quantize(colors=64);frames=[f.quantize(palette=palette,dither=Image.Dither.NONE) for f in frames]
    durations=[110]*len(times);durations[int(np.flatnonzero(times==m.tc)[0])]=800;durations[-1]=1800
    frames[0].save(args.output/'finite_neck.gif',save_all=True,append_images=frames[1:],duration=durations,loop=0,optimize=True)
    snaps=[m.tc-.7,m.tc+.05,m.tc+.7]
    sheet=Image.new('RGB',(1440,360),BG)
    for i,t in enumerate(snaps):sheet.paste(draw_frame(m,t,curves).resize((480,360)),(480*i,0))
    sheet.save(args.output/'sequence.png',optimize=True);draw_frame(m,m.tc+.05,curves).save(args.output/'closeup.png',optimize=True)
    (args.output/'contact.json').write_text(json.dumps(dict(status='KINEMATIC_NOT_DYNAMICAL',parameters=m.par.__dict__,contact_time=m.tc,contact_position=m.xc,time_bracket=m.bracket,contact_gap=m.gap(m.tc),no_wave_contact=models[2].tc,weak_wave_contact=models[1].tc),indent=2,sort_keys=True)+'\n')
    files=[args.output/x for x in ['finite_neck.gif','sequence.png','closeup.png','contact.json']]
    (args.output/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2,sort_keys=True)+'\n')
    print({p.name:p.stat().st_size for p in files})

if __name__=='__main__':main()
