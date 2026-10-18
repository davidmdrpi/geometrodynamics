#!/usr/bin/env python
"""Render option C's smooth density directly, with no surface switch or dissolve."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from neck_wave_contact import ContactNeck, Parameters, H, END, LEVEL

BG = (12,22,34); FG = (222,233,241); MUTED = (150,174,190)
TEAL = (85,215,204); GOLD = (255,195,111)
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
STOPS = np.array([[12,22,34],[46,50,88],[143,66,122],[236,119,143],[255,219,155]])
PALETTE = np.array([np.interp(np.linspace(0,1,256),np.linspace(0,1,5),STOPS[:,i]) for i in range(3)]).T.astype('uint8')


def contour_segments(xs, ys, z, level=LEVEL):
    """Marching squares; resolve saddle cells by the cell-center value."""
    bits = ((z[:-1,:-1]>=level).astype(int)+2*(z[:-1,1:]>=level)+
            4*(z[1:,1:]>=level)+8*(z[1:,:-1]>=level))
    table = {1:[(3,0)],2:[(0,1)],3:[(3,1)],4:[(1,2)],6:[(0,2)],7:[(3,2)],
             8:[(2,3)],9:[(0,2)],11:[(1,2)],12:[(1,3)],13:[(0,1)],14:[(3,0)]}
    result = []
    for j,i in zip(*np.nonzero((bits!=0)&(bits!=15))):
        corners=[(xs[i],ys[j]),(xs[i+1],ys[j]),(xs[i+1],ys[j+1]),(xs[i],ys[j+1])]
        v=[z[j,i],z[j,i+1],z[j+1,i+1],z[j+1,i]]; mask=bits[j,i]
        if mask in (5,10):
            center=sum(v)/4>=level
            pairs=[(0,1),(2,3)] if (mask==5)==center else [(3,0),(1,2)]
        else: pairs=table[mask]
        def edge(k):
            l=(k+1)%4; u=(level-v[k])/(v[l]-v[k])
            return tuple((1-u)*a+u*b for a,b in zip(corners[k],corners[l]))
        result.extend((edge(a),edge(b)) for a,b in pairs)
    return result


def panel(im, model, t, box, extent):
    left,top,width,height=box; xmin,xmax,ymin,ymax=extent
    xs=np.linspace(xmin,xmax,width); ys=np.linspace(ymax,ymin,height)
    z=model.density(xs[None,:],ys[:,None],t)
    im.paste(Image.fromarray(PALETTE[np.rint(z*255).clip(0,255).astype(int)]),(left,top))
    d=ImageDraw.Draw(im)
    for a,b in contour_segments(np.arange(width),np.arange(height),z):
        d.line([(left+a[0],top+a[1]),(left+b[0],top+b[1])],fill=FG,width=1)
    axis=top+height*ymax/(ymax-ymin)
    for x in range(left,left+width,12):d.line([(x,axis),(min(x+5,left+width),axis)],fill=(99,132,144))
    d.rectangle((left,top,left+width,top+height),outline=(54,80,99),width=1)


def draw_frame(model,t,curves):
    im=Image.new('RGB',(960,740),BG);d=ImageDraw.Draw(im)
    font=lambda n,b=False:ImageFont.truetype(BOLD if b else FONT,n)
    d.text((25,18),'Wave contact as a smooth diffuse field',font=font(26,True),fill=FG)
    d.text((25,54),'OPTION C  /  same packet, no imposed split, no cross-fade',font=font(14),fill=GOLD)
    title,note=model.phase(t)
    d.text((25,86),title,font=font(21,True),fill=TEAL)
    d.text((25,116),note,font=font(13),fill=MUTED)
    d.text((930,89),f't = {t:.2f}',font=font(18),fill=FG,anchor='ra')
    d.text((25,147),'FULL NECK: axial section (different axis scales)',font=font(12,True),fill=MUTED)
    d.text((625,147),'FOCUS: equal-scale axial section',font=font(12,True),fill=MUTED)
    panel(im,model,t,(25,172,560,320),(-H,H,-1.08,1.08))
    panel(im,model,t,(625,172,310,310),(-.23,.23,-.23,.23))
    d=ImageDraw.Draw(im)
    for x,val in [(25,'−0.6'),(305,'x = 0'),(585,'+0.6')]:d.text((x,499),val,font=font(12),fill=MUTED,anchor='ma' if x==305 else ('ra' if x==585 else 'la'))
    d.text((625,490),'Dashed: axis ρ = 0. White: ψ = ½.',font=font(12),fill=FG)
    d.text((625,511),'Rotate this section about the axis for 3D.',font=font(12),fill=MUTED)
    # Actual field values, not opacity interpolation between geometries.
    bar=np.repeat(PALETTE[None,:,:],12,axis=0)
    im.paste(Image.fromarray(bar),(25,537));d=ImageDraw.Draw(im)
    d.text((25,553),'ψ = 0',font=font(12),fill=MUTED);d.text((281,553),'1',font=font(12),fill=MUTED,anchor='ra')
    d.text((304,535),'Color is shell density. The field remains smooth;',font=font(13),fill=FG)
    d.text((304,555),'a selected contour can still change topology at a critical point.',font=font(13),fill=MUTED)
    # max on the axis measures filling; after the wave passes it returns to zero.
    px=lambda q:240+q/END*670;py=lambda p:663-p*62
    d.text((25,592),'PEAK AXIS DENSITY',font=font(12,True),fill=MUTED)
    d.text((25,615),'Teal: default packet',font=font(12),fill=TEAL)
    d.text((25,637),'Gray: weak / wave off',font=font(12),fill=MUTED)
    for y,lab in [(0,'0'),(.5,'½'),(1,'1')]:
        d.line([(240,py(y)),(910,py(y))],fill=(44,65,81),width=1)
        d.text((228,py(y)-7),lab,font=font(11),fill=MUTED,anchor='ra')
    for curve,col in zip(curves,[TEAL,(153,177,198),(77,98,118)]):d.line([(px(q),py(p)) for q,p in curve],fill=col,width=2)
    for q in model.field_events:
        d.line([(px(q),594),(px(q),667)],fill=GOLD,width=1)
    d.line([(px(t),596),(px(t),668)],fill=FG,width=2)
    for q in [0,3,6,9]:d.text((px(q),669),str(q),font=font(11),fill=MUTED,anchor='ma')
    d.text((25,709),'Prescribed field, not solved fluid or spacetime dynamics. No mass, energy, or momentum conservation claim.',font=font(12),fill=MUTED)
    return im


def audit(model):
    # A sampled spatial Hessian check, not a certified global maximum.
    xx=np.linspace(-H,H,241)[None,:];rr=np.linspace(0,1.08,217)[:,None]
    maximum=0.
    for t in np.linspace(0,END,91):
        a,b,c,d=model.hessian_components(xx,rr,t)
        maximum=max(maximum,float(np.max(np.sqrt(a*a+b*b+c*c+2*d*d))))
    record=model.record()
    record.update(no_wave_band_crossings=ContactNeck(Parameters(amplitude=0)).band_events,
                  weak_wave_band_crossings=ContactNeck(Parameters(amplitude=.08)).band_events,
                  no_wave_field_crossings=ContactNeck(Parameters(amplitude=0)).field_events,
                  weak_wave_field_crossings=ContactNeck(Parameters(amplitude=.08)).field_events,
                  sampled_hessian_frobenius_max=maximum,
                  sampled_hessian_scaled_max=maximum*model.par.thickness**2,
                  hessian_grid=dict(x_points=241,rho_points=217,time_points=91,x_range=[-H,H],rho_range=[0,1.08],time_range=[0,END]))
    return record


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=Path('docs/images/finite_neck'));args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    m=ContactNeck();models=[m,ContactNeck(Parameters(amplitude=.08)),ContactNeck(Parameters(amplitude=0))]
    curves=[[(t,mm.axis_peak(t)) for t in np.linspace(0,END,181)] for mm in models]
    times=np.unique(np.r_[np.linspace(0,END,61),m.field_events,m.tc])
    frames=[draw_frame(m,float(t),curves).resize((720,555)) for t in times]
    palette=Image.new('RGB',(720,555*6))
    for i,j in enumerate(np.linspace(0,len(times)-1,6).astype(int)):palette.paste(frames[j],(0,555*i))
    palette=palette.quantize(colors=96);frames=[f.quantize(palette=palette,dither=Image.Dither.NONE) for f in frames]
    durations=[120]*len(times);durations[-1]=1600
    frames[0].save(args.output/'finite_neck.gif',save_all=True,append_images=frames[1:],duration=durations,loop=0,optimize=True)
    snaps=[3.2,4.25,5.3];sheet=Image.new('RGB',(1440,370),BG)
    for i,t in enumerate(snaps):sheet.paste(draw_frame(m,t,curves).resize((480,370)),(480*i,0))
    sheet.save(args.output/'sequence.png',optimize=True)
    draw_frame(m,m.field_events[0]+.16,curves).save(args.output/'closeup.png',optimize=True)
    (args.output/'contact.json').write_text(json.dumps(audit(m),indent=2,sort_keys=True)+'\n')
    files=[args.output/x for x in ['finite_neck.gif','sequence.png','closeup.png','contact.json']]
    (args.output/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2,sort_keys=True)+'\n')
    print({p.name:p.stat().st_size for p in files});print(json.dumps(m.record(),indent=2))

if __name__=='__main__':main()
