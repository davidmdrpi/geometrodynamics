#!/usr/bin/env python
"""Klein immersion guide + scalar sphere deformation + diffuse intrinsic throat."""
import argparse,base64,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from geometrodynamics.viz.klein_handles import build_surface,HandleSpec
from klein_release_model import KleinRelease
from render_finite_neck import PALETTE,contour_segments

BG=(12,22,34);FG=(222,233,241);MUTED=(150,174,190);TEAL=(85,215,204);GOLD=(255,195,111);PINK=(227,133,182)
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf';BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
font=lambda n,b=False:ImageFont.truetype(BOLD if b else FONT,n)


class Renderer:
    def __init__(self,m):
        self.m=m
        self.surf=build_surface(handles=[HandleSpec((1,0,0),(0,-1,0),True)],n_sphere=700,n_ring=24,ds_target=.018)
        self.ids=self.surf.tube_ids[0];self.guide=self.surf.points[self.ids]
        self.coord=-self.surf.s[self.ids[:,0]]/self.surf.neck.S*.6
        theta=np.linspace(.001,np.pi-.001,41);phi=np.linspace(0,2*np.pi,81)[:-1]
        tt,pp=np.meshgrid(theta,phi,indexing='ij')
        self.unit=np.stack([np.cos(tt),np.sin(tt)*np.cos(pp),np.sin(tt)*np.sin(pp)],-1).reshape(-1,3)
        self.Y=m.basis(self.unit[:,0]);self.tris=[]
        for i in range(40):
            for j in range(80):
                a=i*80+j;b=i*80+(j+1)%80;c=(i+1)*80+j;d=(i+1)*80+(j+1)%80
                self.tris.extend([(a,b,c),(b,d,c)])
        self.tris=np.array(self.tris)
        self.gtris=[]
        for i in range(len(self.ids)-1):
            for j in range(24):
                a=i*24+j;b=i*24+(j+1)%24;c=(i+1)*24+j;d=(i+1)*24+(j+1)%24
                self.gtris.extend([(a,b,c),(b,d,c)])
        self.gtris=np.array(self.gtris)
        eye=np.array([3.,-5.,2.7]);self.eye=eye/np.linalg.norm(eye)
        self.right=np.cross([0,0,1],self.eye);self.right/=np.linalg.norm(self.right)
        self.up=np.cross(self.eye,self.right)
        self.mu=np.cos(np.linspace(0,np.pi,241));self.lineY=m.basis(self.mu)

    def project(self,p):
        return np.column_stack([325+175*np.asarray(p)@self.right,345-175*np.asarray(p)@self.up])

    def scene(self,im,t):
        m=self.m;q=m.progress(t);a,p=m.state(t);wave=self.Y@a
        # Radial height is a sphere-only analogy; no normal displacement is
        # defined globally on the non-orientable handle.
        theta=np.arccos(np.clip(abs(self.unit[:,0]),0,1))
        indent=-.08*(1-q)*np.exp(-(theta/.3)**4)
        sphere=self.unit*(1+indent+.22*wave[:,None].ravel()/m.peak)[:,None]
        guide=self.guide.copy();center=guide.mean(axis=1,keepdims=True)
        signal=m.transmitted(self.coord,m.throat_clock(t))
        guide=center+(guide-center)*(1-.4*signal[:,None,None]*(1-(self.coord/.6)**2)[:,None,None]**2)
        targets=np.where((self.coord<0)[:,None,None],self.guide[0][None,:,:],self.guide[-1][None,:,:])
        guide=(1-q)*guide+q*targets
        alpha=1-smooth_local((1.25*q-(self.coord/.6)**2)/.25)
        pts=np.vstack([sphere,guide.reshape(-1,3)])
        tri=np.vstack([self.tris,self.gtris+len(sphere)])
        signed=np.r_[wave/m.peak,np.repeat(signal,len(self.ids[0]))]
        faceval=signed[tri].mean(axis=1)
        facealpha=np.r_[np.ones(len(self.tris)),np.repeat(alpha,24)[self.gtris].mean(axis=1)]
        vertices=pts[tri];normal=np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-12)
        light=np.array([-.2,-.6,1.]);light/=np.linalg.norm(light)
        shade=.48+.50*np.abs(normal@light)
        colors=np.empty((len(tri),3));strength=np.minimum(abs(faceval)*2.5,1)
        for i,v in enumerate(faceval):
            base=np.array((53,88,116) if i<len(self.tris) else (176,139,79))
            tint=np.array(TEAL if v>=0 else PINK);colors[i]=(base*(1-strength[i])+tint*strength[i])*shade[i]
        projected=self.project(pts);depth=(vertices@self.eye).mean(axis=1)
        d=ImageDraw.Draw(im,'RGBA')
        for k in np.argsort(depth):
            if facealpha[k]<.002:continue
            d.polygon([tuple(v) for v in projected[tri[k]]],fill=(*np.clip(colors[k],0,255).astype(int),int(facealpha[k]*255)))
        d=ImageDraw.Draw(im)
        for sign,label in [(1,'A'),(-1,'B / antipode')]:
            ph=np.linspace(0,2*np.pi,81)
            ring=np.column_stack([np.full(len(ph),sign*np.cos(.3)),np.sin(.3)*np.cos(ph),np.sin(.3)*np.sin(ph)])
            if q<.999:d.line([tuple(v) for v in self.project(ring)],fill=GOLD,width=2)
            pos=self.project(np.array([[sign*1.24,0,0]]))[0]
            d.text(tuple(pos),label,font=font(14,True),fill=GOLD,anchor='mm')
        d.text((28,572),'Curved Klein handle: identification guide, not physical length',font=font(13),fill=GOLD)
        d.text((28,592),'Sphere height and color show the computed signed scalar field.',font=font(13),fill=MUTED)

    def frame(self,t):
        m=self.m;im=Image.new('RGB',(1100,800),BG);d=ImageDraw.Draw(im)
        d.text((24,17),'Sphere wave → throat contact → release → antipode',font=font(25,True),fill=FG)
        d.text((24,54),'KLEIN IDENTIFICATION  /  intrinsic throat shown separately  /  irreversible release is a chosen rule',font=font(13),fill=GOLD)
        q=m.progress(t)
        if t<m.focus-.12:phase='Incoming deformation converges on mouth A'
        elif t<m.tc:phase='Mouth signal propagates down the intrinsic throat'
        elif q<1:phase='Diffuse contact; the two throat halves withdraw'
        elif t<m.tc+m.pulse_duration+2.5:phase='Released wave travels over the sphere'
        elif t<m.tc+m.pulse_duration+3.7:phase='Released packet refocuses near the antipodal mouths'
        else:phase='Standing-mode superposition continues its return cycle'
        d.text((24,87),phase,font=font(20,True),fill=TEAL)
        d.text((1074,87),f't = {t:.2f}',font=font(18),fill=FG,anchor='ra')
        self.scene(im,t);d=ImageDraw.Draw(im)
        d.text((738,135),'INTRINSIC THROAT',font=font(15,True),fill=FG)
        d.text((738,157),'Smooth density ψ; white contour = ½',font=font(12),fill=MUTED)
        w,h=330,280;extent=.3+.35*q;vertical=.25+.8*q
        xs=np.linspace(-extent,extent,w);ys=np.linspace(vertical,-vertical,h)
        values=m.density(xs[None,:],ys[:,None],t)
        im.paste(Image.fromarray(PALETTE[np.rint(values*255).clip(0,255).astype(int)]),(738,182));d=ImageDraw.Draw(im)
        for a,b in contour_segments(np.arange(w),np.arange(h),values):d.line([(738+a[0],182+a[1]),(738+b[0],182+b[1])],fill=FG,width=1)
        d.rectangle((738,182,1068,462),outline=(54,80,99))
        d.text((738,474),'Contact depends on the incoming mouth signal.',font=font(12),fill=MUTED)
        d.text((738,494),'Withdrawal after contact is prescribed.',font=font(12),fill=GOLD)
        if q>.999:d.text((903,322),'Connection released',font=font(17,True),fill=TEAL,anchor='mm')
        ledger=m.ledger(t)
        d.text((738,532),'ASSIGNED RELEASE BUDGET',font=font(13,True),fill=MUTED)
        d.text((738,557),f'Reservoir remaining   {ledger["reservoir"]:+.5f}',font=font(15),fill=GOLD)
        d.text((738,581),f'Added field energy     {ledger["work"]:+.5f}',font=font(15),fill=TEAL)
        d.text((738,605),'Includes interference; not derived neck energy.',font=font(11),fill=MUTED)
        # Meridian traces make the emitted component visible within the total.
        d.line([(25,632),(1075,632)],fill=(45,66,82))
        d.text((25,644),'SPHERE MERIDIAN',font=font(12,True),fill=MUTED)
        d.text((25,666),'White: total field',font=font(12),fill=FG)
        d.text((25,686),'Teal: emitted component',font=font(12),fill=TEAL)
        a,p=m.state(t);ae,pe=m.state(t,emitted_only=True)
        px=np.linspace(240,1068,241)
        for data,col in [(self.lineY@a,FG),(self.lineY@ae,TEAL)]:
            d.line(list(zip(px,701-34*data/m.peak)),fill=col,width=2)
        d.line([(240,701),(1068,701)],fill=(54,80,99))
        for x,lab in [(240,'A (θ = 0)'),(654,'equator'),(1068,'B (θ = π)')]:d.text((x,739),lab,font=font(12),fill=MUTED,anchor='ma' if x==654 else ('ra' if x==1068 else 'la'))
        d.text((25,776),'Illustrative mechanism: solved sphere modes + delayed throat readout + imposed release. No coupled GR/fluid claim.',font=font(12),fill=MUTED)
        return im


def smooth_local(q):
    q=np.clip(q,0,1);return q*q*q*(10+q*(-15+6*q))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('docs/images/klein_release'));args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    m=KleinRelease();r=Renderer(m)
    times=np.unique(np.r_[np.linspace(0,m.end,65),m.focus,m.tc,m.tc+.08,m.tc+.16,m.tc+.24,m.tc+.32,m.tc+.4,m.tc+.5,m.tc+.6,m.tc+m.withdrawal,m.tc+m.pulse_duration+np.pi])
    frames=[r.frame(float(t)).resize((660,480)) for t in times]
    pal=Image.new('RGB',(660,480*6))
    for i,j in enumerate(np.linspace(0,len(times)-1,6).astype(int)):pal.paste(frames[j],(0,i*480))
    pal=pal.quantize(colors=64);indexed=[f.quantize(palette=pal,dither=Image.Dither.NONE) for f in frames]
    durations=[120]*len(times);durations[-1]=1600
    for event in [m.tc,m.tc+.4,m.tc+m.pulse_duration+np.pi]:durations[np.argmin(abs(times-event))]=700
    indexed[0].save(args.output/'klein_release.gif',save_all=True,append_images=indexed[1:],duration=durations,loop=0,optimize=True)
    for name,t in [('contact',m.tc+.12),('retraction',m.tc+.4),('antipode',m.tc+m.pulse_duration+np.pi)]:r.frame(t).save(args.output/(name+'.png'),optimize=True)
    (args.output/'report.json').write_text(json.dumps(m.report(),indent=2,sort_keys=True)+'\n')
    # Browser-free, self-contained scrubber using the exact rendered frames.
    images=[]
    import io
    stage_times=[0.,.85,m.focus,m.tc+.12,m.tc+.32,m.tc+m.withdrawal,m.tc+m.pulse_duration/2+np.pi/2,m.tc+m.pulse_duration/2+np.pi,m.end]
    for t in stage_times:
        buf=io.BytesIO();r.frame(t).save(buf,format='WEBP',quality=65,method=6);images.append('data:image/webp;base64,'+base64.b64encode(buf.getvalue()).decode())
    html='''<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Klein throat release</title><style>body{max-width:1000px;margin:24px auto;padding:12px;background:#0c1622;color:#dee9f1;font:16px system-ui}img,input{width:100%}button{padding:10px;margin:8px;background:#214655;color:white;border:1px solid #55d7cc}p{line-height:1.5;color:#a6bdce}</style><h1>Klein throat release — nine-stage viewer</h1><p>The curved handle illustrates an extradimensional identification. The inset shows the intrinsic diffuse field. Separation and withdrawal are imposed; the sphere modes and assigned pulse-work budget are computed.</p><img id="frame" alt="Klein-handle collapse and sphere wave"><button id="play">Play</button><label>Frame <output id="out"></output><input id="scrub" type="range" min="0" value="0"></label><p>White: total sphere deformation. Teal: emitted component. The traveling packet is a superposition of standing modes; it is not a stationary single mode.</p><script>const frames=FRAMES,times=TIMES;let i=0,playing=false;const img=document.getElementById('frame'),scrub=document.getElementById('scrub'),out=document.getElementById('out'),play=document.getElementById('play');scrub.max=frames.length-1;function draw(){img.src=frames[i];scrub.value=i;out.textContent='t = '+times[i].toFixed(3)}scrub.oninput=()=>{i=+scrub.value;playing=false;play.textContent='Play';draw()};play.onclick=()=>{playing=!playing;play.textContent=playing?'Pause':'Play'};setInterval(()=>{if(playing){i=(i+1)%frames.length;draw()}},900);draw();</script>'''.replace('FRAMES',json.dumps(images)).replace('TIMES',json.dumps(stage_times))
    viewer=args.output.parent.parent/'visualizations'/'klein_release.html';viewer.write_text(html)
    files=[p for p in args.output.iterdir() if p.name!='sha256.json']
    (args.output/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2,sort_keys=True)+'\n')
    print(json.dumps(m.report(),indent=2));print({p.name:p.stat().st_size for p in files});print('viewer bytes',viewer.stat().st_size)

if __name__=='__main__':main()
