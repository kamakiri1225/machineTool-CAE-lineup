"""Plot solver FRFs and animate complex forced response (not eigenmodes)."""
from pathlib import Path
import sys,re,json
import numpy as np
from PIL import Image,ImageDraw,ImageFont
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'tools'))
from plot_modes import plt
import mpc_springs as ms
import vtu_tools as vt
import pyvista as pv
from frequency_response import SELECT
R=D/'results/frf'; R.mkdir(parents=True,exist_ok=True)

def ccx_dat(path):
    lines=path.read_text().splitlines(); rec=[]
    for i,s in enumerate(lines):
        m=re.match(r'\s*displacements .* and time\s+([\d.E+-]+)',s)
        if m: rec.append((float(m[1]),np.array(list(map(float,lines[i+2].split()[1:4])))))
    steps=[]; cur=[]
    for i in range(0,len(rec),2):
        f,u=rec[i]; z=u+1j*rec[i+1][1]
        if cur and f<cur[-1][0]: steps.append(cur);cur=[]
        cur.append((f,z))
    steps.append(cur)
    assert len(steps)==3, len(steps)
    return [(np.array([x[0] for x in s]),np.array([x[1] for x in s])) for s in steps]

def read_front(path):
    lines=path.read_text().splitlines(); k=lines.index('*data'); nn=int(lines[k+1].split()[0]); nc=int(lines[k+2].split()[0]); nd=sum(map(int,lines[k+3].split())); labels=lines[k+4:k+4+nc]; offsets=np.cumsum([0]+list(map(int,lines[k+3].split()))); ir=next(j for j,l in enumerate(labels) if l.lower()=='displacement_real'); ii=next(j for j,l in enumerate(labels) if l.lower()=='displacement_imag'); k=k+4+nc
    result={}
    for _ in range(nn):
        n=int(lines[k]);k+=1; values=[]
        while len(values)<nd: values.extend(map(float,lines[k].split())); k+=1
        result[n]=np.array(values[offsets[ir]:offsets[ir]+3])+1j*np.array(values[offsets[ii]:offsets[ii]+3])
    return result

def selected_frd(path,targets):
    # Streaming fixed-width reader: only retain 3 complex responses, discard eigenmodes and other sweep points.
    chosen=[None]*3; best=[float('inf')]*3; sweep=-1; prev=-1; real=None
    with path.open() as stream:
        for line in stream:
            if not line.startswith('  100CL'): continue
            freq=float(line[12:24]); name=next(stream)[5:13].strip()
            if name not in ('DISP','DISPI'): continue
            data={}
            for row in stream:
                if row.startswith(' -1'): break
            while row.startswith(' -1'):
                data[int(row[3:13])]=np.array([float(row[13:25]),float(row[25:37]),float(row[37:49])]); row=next(stream)
            # eigenmode frequencies increase, then first sweep returns to 50Hz
            if freq<prev-1e-4: sweep+=1
            prev=freq
            if sweep<0 or sweep>2: continue
            if real is None: real=(freq,data)
            else:
                assert abs(real[0]-freq)<1e-3
                err=abs(freq-targets[sweep])
                if err<best[sweep]:
                    chosen[sweep]=(freq,{n:real[1][n]+1j*v for n,v in data.items()});best[sweep]=err
                real=None
    assert all(x is not None for x in chosen),(best,sweep)
    return chosen

def animate(kind,dof,model,cc,front):
    grid,ids=vt.make_grid(type('M',(),{'parts':model.parts,'springs':[]})())
    grid.point_data['index']=np.arange(len(ids))
    surf=grid.extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface')
    ind=surf.point_data['index'].astype(int); original=surf.points.copy()
    uc=np.array([cc[1][n] for n in ids])[ind]; uf=np.array([front[n] for n in ids])[ind]
    maxu=max(np.linalg.norm(abs(uc),axis=1).max(),np.linalg.norm(abs(uf),axis=1).max()); scale=65/maxu
    font=ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),18)
    frames=[]; p=pv.Plotter(shape=(1,2),off_screen=True,window_size=(1250,650)); surfaces=[]
    for panel,(name,u) in enumerate([('CalculiX',uc),('FrontISTR',uf)]):
        p.subplot(0,panel);p.set_background('white')
        p.add_mesh(surf,style='wireframe',opacity=.12,color='gray')
        s=surf.copy();s.point_data['amplitude']=np.linalg.norm(abs(u),axis=1)*1000
        p.add_mesh(s,scalars='amplitude',clim=(0,maxu*1000),cmap='turbo',show_edges=False,show_scalar_bar=False)
        p.camera_position=[(2100,-1700,1400),(500,380,400),(0,0,1)]
        p.camera.zoom(0.80); p.add_axes();surfaces.append(s)
    for phase in np.linspace(0,2*np.pi,32,endpoint=False):
        for s,u in zip(surfaces,[uc,uf]): s.points=original+scale*np.real(u*np.exp(1j*phase))
        p.render()
        img=Image.fromarray(p.screenshot(return_img=True)); draw=ImageDraw.Draw(img)
        draw.rectangle((0,0,1250,70),fill='white')
        draw.text((12,5),f'{kind.upper()}・{"XYZ"[dof-1]}方向 1 N 調和加振  位相 {phase*180/np.pi:.0f}°   変形倍率 {scale:.0f}倍',font=font,fill='black')
        draw.text((12,34),f'CalculiX {cc[0]:.3f} Hz',font=font,fill='black')
        draw.text((640,34),f'FrontISTR {SELECT[kind][dof-1]:.3f} Hz',font=font,fill='black')
        draw.rectangle((0,610,1250,650),fill='white')
        draw.text((20,616),'色：変位振幅 [µm / 1 N]',font=font,fill='black')
        for x in range(240):
            color=tuple(int(255*v) for v in plt.get_cmap('turbo')(x/239)[:3])
            draw.line((300+x,622,300+x,638),fill=color)
        draw.text((270,613),'0',font=font,fill='black')
        draw.text((550,613),f'{maxu*1000:.4f}',font=font,fill='black')
        draw.text((770,616),'灰色線：変形前の形状',font=font,fill='black')
        frames.append(img)
    p.close(); frames[0].save(R/f'{kind}_{dof}.gif',save_all=True,append_images=frames[1:],duration=65,loop=0)
    frames[0].save(R/f'{kind}_{dof}.png')

def main():
    rows=['| 接続 | 方向 | CalculiXピーク Hz | 振幅 µm/N | FrontISTRピーク Hz | 振幅 µm/N |','|---|---|---:|---:|---:|---:|']; checks=[]
    for kind in SELECT:
        base=D/'frf'/kind; cc=ccx_dat(base/'ccx/response.dat')
        shapes=None if '--plots-only' in sys.argv else selected_frd(base/'ccx/response.frd',SELECT[kind]);model=ms.build(kind)
        fig,axs=plt.subplots(3,1,figsize=(11,10),sharex=True)
        for i,(f,u) in enumerate(cc):
            dof=i+1; cf=np.linalg.norm(abs(u),axis=1)*1000
            data=[]
            for ln in (base/f'fistr/sweep_{dof}.log').read_text().splitlines():
                m=re.match(r'\s*([\d.E+-]+)\s+\[Hz\]\s*:\s*([\d.E+-]+)\s*$',ln)
                if m:data.append(tuple(map(float,m.groups())))
            ff,fa=np.array(data).T;fa*=1000
            ax=axs[i];ax.plot(f,cf,label='CalculiX',lw=2);ax.plot(ff,fa,'--',label='FrontISTR',lw=1.5)
            ax.set_title(f'{kind.upper()}：テーブル上面中央に {"XYZ"[i]}方向 1 N → 同じ点の3成分振幅ノルム',loc='left');ax.set_ylabel('変位振幅 [µm/N]');ax.grid(alpha=.3,which='both');ax.legend()
            j=cf.argmax();k=fa.argmax();rows.append(f'| {kind.upper()} | {"XYZ"[i]} | {f[j]:.3f} | {cf[j]:.5f} | {ff[k]:.3f} | {fa[k]:.5f} |')
            np.savetxt(R/f'{kind}_{dof}_ccx.csv',np.column_stack([f,u.real,u.imag,cf]),delimiter=',',header='frequency_Hz,Ux_real_mm,Uy_real_mm,Uz_real_mm,Ux_imag_mm,Uy_imag_mm,Uz_imag_mm,norm_um_per_N')
            np.savetxt(R/f'{kind}_{dof}_fistr.csv',np.column_stack([ff,fa]),delimiter=',',header='frequency_Hz,norm_um_per_N')
            ax.set_ylim(bottom=0)
            if '--plots-only' in sys.argv: continue
            front=read_front(base/f'fistr/shape_{dof}_dyna.res.0.1')
            # Confirm the shape result and the selected-frequency scalar log use the same solution.
            monitor=np.linalg.norm(abs(front[model.monitor]))
            log=(base/f'fistr/shape_{dof}.log').read_text();scalar=float(re.search(r'\[Hz\]\s*:\s*([\d.E+-]+)\s*\n',log)[1])
            assert np.isclose(monitor,scalar,rtol=1e-5),(monitor,scalar)
            # Match FRD vector to corresponding DAT output at the selected CCX frequency.
            q=np.argmin(abs(f-shapes[i][0])); err=np.linalg.norm(shapes[i][1][model.monitor]-u[q])/max(np.linalg.norm(u[q]),1e-30)
            assert err<2e-3,err
            checks.append({'kind':kind,'dof':dof,'ccx_gif_frequency':shapes[i][0],'front_gif_frequency':SELECT[kind][i],'ccx_frd_dat_relative_error':err,'front_shape_log_relative_error':abs(monitor/scalar-1)})
            animate(kind,dof,model,shapes[i],front)
        axs[-1].set_xlabel('加振周波数 [Hz]');fig.tight_layout();fig.savefig(R/f'{kind}_frf.png');plt.close(fig)
    (R/'summary.md').write_text('\n'.join(rows)+'\n')
    if checks: (R/'checks.json').write_text(json.dumps(checks,indent=2))
    print('\n'.join(rows))

if __name__ == '__main__':
    main()
