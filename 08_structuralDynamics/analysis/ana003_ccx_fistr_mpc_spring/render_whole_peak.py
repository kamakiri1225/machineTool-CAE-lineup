"""Assembled whole-machine forced vibration above the matching FRF peak."""
from pathlib import Path
import sys,io
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont,ImageChops
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'tools'))
import mpc_springs as ms
import vtu_tools as vt
from plot_modes import plt
from post_frequency import read_front
foundation='--foundation' in sys.argv
R=D/'foundation/results' if foundation else D/'results/frf'
source=D/'foundation/fistr' if foundation else D/'frf/rbe2/fistr'
font=lambda n:ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),n)
if foundation:
 from foundation_model import build
 m=build()
else: m=ms.build('rbe2')
g,ids=vt.make_grid(type('M',(),{'parts':m.parts,'springs':[]})());g.point_data['node_id']=ids
import json
selected=json.loads((R/'selected.json').read_text()) if foundation else [183.,454.,625.]
for dof,freq in enumerate(selected,1):
 data=read_front(source/f'shape_{dof}_dyna.res.0.1')
 u=np.array([data[n] for n in ids]);amplitude=np.linalg.norm(abs(u),axis=1);scale=85/amplitude.max()
 grid=g.copy();grid.point_data['U_real']=u.real;grid.point_data['U_imag']=u.imag;grid.point_data['amplitude_um']=amplitude*1000;grid.save(R/f'whole_response_{dof}.vtu')
 curve=np.loadtxt(R/f'rbe2_{dof}_fistr_components.csv',delimiter=',');amps=abs(data[m.monitor])*1000
 direction='XYZ'[dof-1];labels=[f"H{axis.lower()}{direction.lower()} : {direction}加振 → {axis}変位" for axis in 'XYZ'];colors=['#2070b4','#d67b13','#248347']
 fig,ax=plt.subplots(figsize=(10,3.2),dpi=120)
 for comp,(label,color) in enumerate(zip(labels,colors)):
  ax.plot(curve[:,0],curve[:,comp+1],label=label,color=color,lw=1.8)
  ax.plot(freq,amps[comp],'o',color=color,ms=7,markeredgecolor='red')
 ax.axvline(freq,color='#d52d34',ls='--')
 ax.annotate(f'上段のアニメーション：{freq:.0f} Hz',xy=(freq,amps.max()),xytext=(.53,.87),textcoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#d52d34'},fontsize=10,color='#b22222')
 ax.set_xlabel('加振周波数 [Hz]');ax.set_ylabel('同じMON点の方向別コンプライアンス [µm/N]');ax.set_xlim(5 if foundation else 50,750)
 ax.set_ylim(0,curve[:,1:].max()*1.18);ax.grid(alpha=.3,which='both');ax.legend(loc='lower left',fontsize=9)
 fig.tight_layout();buf=io.BytesIO();fig.savefig(buf,format='png');plt.close(fig);buf.seek(0);chart=Image.open(buf).convert('RGB').resize((1200,380))
 p=pv.Plotter(off_screen=True,window_size=(1200,720));p.set_background('white');animated=[]
 for part,color in [(0,'#a7b3c2'),(1,'#328ccb')]:
  surf=grid.extract_cells(np.where(grid.cell_data['part']==part)[0]).extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface')
  base=surf.points.copy();du=surf.point_data['U_real']+1j*surf.point_data['U_imag']
  p.add_mesh(surf.copy(),style='wireframe',color='#737373',opacity=.18,line_width=.4)
  moving=surf.copy();p.add_mesh(moving,color=color,show_edges=False,smooth_shading=False)
  animated.append((moving,base,du))
 if foundation: p.add_mesh(pv.Plane(center=(500,450,-80),direction=(0,0,1),i_size=1150,j_size=1050),color='#d3dfcc',opacity=.75)
 monitor0=np.array(next(c[1:] for part in m.parts for c in part.coords if c[0]==m.monitor),dtype=float)
 force_dir=np.eye(3)[dof-1]
 marker=pv.PolyData(monitor0[None,:]);p.add_mesh(marker,color='#f01824',point_size=19,render_points_as_spheres=True)
 arrow=pv.Arrow(start=monitor0-force_dir*200,direction=force_dir,scale=185,tip_length=.23,tip_radius=.08,shaft_radius=.025)
 arrowbase=arrow.points.copy();p.add_mesh(arrow,color='#f01824',lighting=False)
 p.camera_position=[(2100,-1800,1500),(500,400,480),(0,0,1)];p.enable_parallel_projection();p.reset_camera();p.camera.zoom(.82);p.add_axes()
 frames=[]
 for phase in np.linspace(np.pi/2,2*np.pi+np.pi/2,40,endpoint=False):
  for mesh,base,du in animated:mesh.points=base+scale*np.real(du*np.exp(1j*phase))
  marker_offset=scale*np.real(data[m.monitor]*np.exp(1j*phase))
  marker.points=(monitor0+marker_offset)[None,:];arrow.points=arrowbase+marker_offset
  p.render()
  top=Image.fromarray(p.screenshot(return_img=True));out=Image.new('RGB',(1200,1100),'white');out.paste(top,(0,0));out.paste(chart,(0,720));draw=ImageDraw.Draw(out)
  draw.rectangle((0,0,1200,83),fill='white');draw.text((15,8),f'{"基礎ばね支持・" if foundation else ""}FrontISTR：テーブル上面中央MONを+{"XYZ"[dof-1]}方向へ1 N加振 / {freq:.0f} Hz',font=font(24),fill='black')
  draw.text((15,46),f"観測：同じMONのX・Y・Z変位（クロス方向を含む3成分）　位相 {phase*180/np.pi%360:.0f}°",font=font(19),fill='black')
  p.renderer.SetWorldPoint(*(monitor0+marker_offset),1.);p.renderer.WorldToDisplay();sx,sy,_=p.renderer.GetDisplayPoint();sy=720-sy
  tx,ty=min(sx+60,840),max(sy-85,110)
  draw.line((sx,sy,tx,ty+27),fill='#d52d34',width=2)
  draw.rectangle((tx-5,ty-3,tx+340,ty+52),fill='white')
  draw.text((tx,ty),'加振点 MON（上面中央）',font=font(20),fill='#c52222')
  draw.text((tx,ty+27),f"+{'XYZ'[dof-1]}方向・振幅1 N／同じ点で変位観測",font=font(16),fill='#c52222')
  draw.rectangle((0,681,1200,720),fill='white');draw.text((15,688),f'形状の切出し・分解なし。見えるよう変位のみ {scale:.0f}倍、1周期をスロー再生。',font=font(18),fill='black')
  frames.append(out)
 p.close()
 path=R/f'whole_peak_response_{dof}.gif';frames[0].save(path,save_all=True,append_images=frames[1:],duration=75,loop=0)
 strip=Image.new('RGB',(1200,720),'white')
 for j,k in enumerate([0,10,20,30]):strip.paste(frames[k].crop((0,0,1200,720)).resize((600,360)),((j%2)*600,(j//2)*360))
 strip.save(R/f'whole_peak_response_{dof}_phases.png')
 assert ImageChops.difference(frames[0].crop((0,90,1200,680)),frames[20].crop((0,90,1200,680))).getbbox()
 print(path.name, '40 frames: whole-part motion verified',flush=True)
