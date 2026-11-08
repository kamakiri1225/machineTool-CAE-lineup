"""Reuse stored 80-mode eigenbasis for selected cutting-force response shapes."""
from pathlib import Path
import sys,shutil,io,numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont,ImageChops
from analyze import read_set
from vmc_model import build,CONTACT
D=Path(__file__).resolve().parent;R=D/'results';case=D/'foundation';sys.path[:0]=[str(D.parent/'tools'),str(D.parent/'ana003_ccx_fistr_mpc_spring')]
from frequency_response import execute,ALPHA,BETA
from frd import read_frd
from plot_modes import plt
import vtu_tools as vt
m=build('foundation');f,w=read_set(case/'relative.dat','WORK')[0];ft,t=read_set(case/'relative.dat','TOOL')[0];q=w-t
base=(case/'relative.inp').read_text().split('*STEP')[0];records=[]
for name,band in [('low',(100,600)),('high',(800,3000))]:
 valid=(f>=band[0])&(f<=band[1]);j=np.where(valid)[0][np.argmax(abs(q[valid,0]))];freq=float(f[j]);job='shape_'+name
 text=base+'\n'.join(['*STEP','*STEADY STATE DYNAMICS',f'{freq}, {freq+.001}, 2, 1.0','*MODAL DAMPING, RAYLEIGH',f',,{ALPHA},{BETA}','*CLOAD, OP=NEW','WORK, 1, 1.0','TOOL, 1, -1.0','*NODE PRINT, NSET=WORK','U','*NODE PRINT, NSET=TOOL','U','*NODE FILE','U','*END STEP'])+'\n'
 (case/f'{job}.inp').write_text(text);shutil.copyfile(case/'relative.eig',case/f'{job}.eig');execute([Path.home()/'opt/ccx.sh','-i',job],case,f'{job}.log')
 coords,blocks=read_frd(case/f'{job}.frd');real=next(b for b in blocks if b[0]=='DISP');imag=next(b for b in blocks if b[0]=='DISPI');data={n:np.array(real[2][n])+1j*np.array(imag[2][n]) for n in real[2]}
 # Compare the selected shape's printed contact response against the sweep result.
 lines=(case/f'{job}.dat').read_text().splitlines();by={'WORK':[],'TOOL':[]}
 for ii,line in enumerate(lines):
  for key in by:
   if f'for set {key} and time' in line:by[key].append(np.array(list(map(float,lines[ii+2].split()[1:4]))))
 sw=by['WORK'][0]+1j*by['WORK'][1];st=by['TOOL'][0]+1j*by['TOOL'][1];error=np.linalg.norm((sw-st)-q[j])/np.linalg.norm(q[j]);assert error<.003,error
 grid,ids=vt.make_grid(type('M',(),{'parts':m.parts,'springs':[]})());grid.point_data['node_id']=ids;grid.cell_data['component']=m.element_parts;u=np.array([data[n] for n in ids]);grid.point_data['U_real']=u.real;grid.point_data['U_imag']=u.imag;grid.save(R/f'cutting_{name}.vtu');scale=65/np.linalg.norm(abs(u),axis=1).max()
 fig,ax=plt.subplots(figsize=(10,3.1),dpi=120)
 for i,color in enumerate(['#246ba5','#cf7a21','#29944a']):
  ax.plot(f,abs(q[:,i])*1000,label=f'G{"xyz"[i]}x：相対{"XYZ"[i]}変位 / X切削力',color=color)
  ax.plot(freq,abs(q[j,i])*1000,'o',color=color,ms=6)
 ax.axvline(freq,color='red',ls='--');ax.set_xlabel('周波数 [Hz]');ax.set_ylabel('相対コンプライアンス [µm/N]');ax.set_ylim(bottom=0);ax.grid(alpha=.3);ax.legend(fontsize=9);fig.tight_layout();buffer=io.BytesIO();fig.savefig(buffer,format='png');plt.close(fig);buffer.seek(0);chart=Image.open(buffer).convert('RGB').resize((1200,370))
 p=pv.Plotter(off_screen=True,window_size=(1200,730));p.set_background('white');meshes=[]
 for i,color in enumerate(['#b5bbc3','#4f94c3','#a8cddc','#c93338','#e3ad46']):
  surf=grid.extract_cells(np.where(np.array(m.element_parts)==i)[0]).extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface');pts=surf.points.copy();du=surf.point_data['U_real']+1j*surf.point_data['U_imag'];p.add_mesh(surf.copy(),style='wireframe',color='gray',opacity=.1);moving=surf.copy();p.add_mesh(moving,color=color);meshes.append((moving,pts,du))
 p.add_mesh(pv.Plane(center=(500,450,-70),direction=(0,0,1),i_size=1150,j_size=1050),color='#d9e4cd',opacity=.6)
 force_meshes=[]
 for force,color in [(1,'#dc2525'),(-1,'#248040')]:
  arrow=pv.Arrow(start=CONTACT,direction=(force,0,0),scale=120,tip_length=.24,tip_radius=.07,shaft_radius=.025);origin=arrow.points.copy();p.add_mesh(arrow,color=color,lighting=False)
  mark=pv.PolyData(CONTACT[None,:]);p.add_mesh(mark,color=color,point_size=13,render_points_as_spheres=True);force_meshes.append((arrow,origin,mark,sw if force==1 else st))
 p.camera_position=[(2100,-1800,1600),(500,450,450),(0,0,1)];p.enable_parallel_projection();p.reset_camera();p.camera.zoom(.85);p.add_axes();font=lambda n:ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),n)
 frames=[]
 for phase in np.linspace(0,2*np.pi,40,endpoint=False):
  for mesh,pts,du in meshes:mesh.points=pts+scale*np.real(du*np.exp(1j*phase))
  for arrow,origin,mark,delta in force_meshes:
   offset=scale*np.real(delta*np.exp(1j*phase));arrow.points=origin+offset;mark.points=(CONTACT+offset)[None,:]
  p.render();top=Image.fromarray(p.screenshot(return_img=True));im=Image.new('RGB',(1200,1100),'white');im.paste(top,(0,0));im.paste(chart,(0,730));dr=ImageDraw.Draw(im);dr.rectangle((0,0,1200,95),fill='white');dr.text((15,8),f'立型MC・基礎ばね支持：X方向の切削力対 / {freq:.1f} Hz',font=font(23),fill='black');dr.text((15,43),'入力：ワーク+1 N／刃先−1 N　観測：同じ切削位置の相対XYZ変位',font=font(19),fill='black');dr.text((15,72),f'位相 {phase*180/np.pi:.0f}°　変位 {scale:.0f}倍　赤：ワーク+X、緑：刃先−X（スロー再生）',font=font(18),fill='black');frames.append(im)
 p.close();frames[0].save(R/f'cutting_{name}.gif',save_all=True,append_images=frames[1:],duration=70,loop=0)
 assert ImageChops.difference(frames[0].crop((180,110,1050,700)),frames[20].crop((180,110,1050,700))).getbbox()
 records.append({'band':name,'frequency_Hz':freq,'shape_sweep_relative_error':float(error),'relative_amplitude_um_N':(abs(q[j])*1000).tolist(),'work_amplitude_um_N':(abs(w[j])*1000).tolist(),'tool_amplitude_um_N':(abs(t[j])*1000).tolist()})
 import json
 (R/'cutting_shapes.json').write_text(json.dumps(records,indent=2));print(records[-1],flush=True)
