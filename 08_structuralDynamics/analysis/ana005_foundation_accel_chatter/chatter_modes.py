"""Boundary feedback shapes, structural-mode contributions, and side-by-side SLD."""
from pathlib import Path
import sys,json,io
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont,ImageChops
D=Path(__file__).resolve().parent;R=D/'results';ROOT=D.parent.parent
sys.path[:0]=[str(ROOT/'skills/chatter-mode-analysis/scripts'),str(D.parent/'tools')]
from boundary_mode import calculate,transfer
from vmc_model import build,CONTACT
from plot_modes import plt
import vtu_tools as vt
with np.load(R/'foundation_modal_basis.npz') as z:basis={k:z[k] for k in z.files}
# Check recovered representative-point modal fields against actual solver outputs.
cols=[np.loadtxt(R/f'foundation_relative_{j}.csv',delimiter=',') for j in (1,2)]
errors=[]
for index in [0,50,150,len(cols[0])//2,len(cols[0])-1]:
 f=cols[0][index,0];G=np.stack([a[index,1:3]+1j*a[index,4:6] for a in cols],axis=1);error=np.linalg.norm(transfer(basis,f)-G)/np.linalg.norm(G);assert error<2e-4,error;errors.append(float(error))
summary=json.loads((R/'milling_checks.json').read_text())['foundation'];states=[]
for speed in [3000,6000,10000]:
 c=summary[f'at_{speed}_rpm'];condition={'rpm':c['evaluated_rpm'],'depth_mm':c['axial_depth_mm'],'chatter_hz':c['chatter_frequency_Hz'],'kt':1000.,'teeth':4,'B0':[[.3,1.],[-1.,.3]]}
 report,U,coeff,share,energy=calculate(basis,condition,refine=True);assert report['loop_relative_residual']<1e-8
 (R/f'chatter_condition_{speed}.json').write_text(json.dumps(condition,indent=2));(R/f'chatter_mode_{speed}.json').write_text(json.dumps(report,indent=2));np.savez_compressed(R/f'chatter_mode_{speed}.npz',U=U,coeff=coeff,share=share,energy=energy)
 np.savetxt(R/f'chatter_mode_{speed}.csv',np.column_stack([np.arange(1,len(coeff)+1),basis['frequency_hz'],coeff.real,coeff.imag,share,energy]),delimiter=',',header='mode,natural_hz,coefficient_real,coefficient_imag,relative_amplitude_share,modal_kinetic_share')
 states.append((speed,report,U))
selected=states[1];speed,report,U=selected;c=report['condition'];top=report['dominant_modes'][:2]
# Plot the same selected condition, distinguishing interpolated SLD and modal refinement.
fig,axs=plt.subplots(2,1,figsize=(9.5,12),dpi=120,sharex=True)
for name,label in [('fixed','底面固定'),('foundation','基礎ばね'),('foundation_table_stiff','基礎ばね・テーブル結合10倍')]:
 a=np.loadtxt(R/f'milling_{name}.csv',delimiter=',');axs[0].plot(a[:,0],a[:,1],label=label);axs[1].plot(a[:,0],a[:,2],label=label)
for ax,y in [(axs[0],c['depth_mm']),(axs[1],c['chatter_hz'])]:ax.plot(c['rpm'],y,'*',color='#c8262d',ms=15);ax.grid(alpha=.3);ax.legend(fontsize=9)
axs[0].set_ylabel('軸方向切込み限界 ap [mm]');axs[0].set_ylim(bottom=0);axs[1].set_ylabel('安定限界でのびびり周波数 fch [Hz]');axs[1].set_xlabel('主軸回転数 [min⁻¹]')
axs[0].annotate(f"右のびびり形状\nn={c['rpm']:.1f} min⁻¹\nap={c['depth_mm']:.3f} mm",xy=(c['rpm'],c['depth_mm']),xytext=(.45,.82),textcoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#c8262d'},color='#a9242b',fontsize=11)
axs[1].annotate(f"fch={c['chatter_hz']:.1f} Hz",xy=(c['rpm'],c['chatter_hz']),xytext=(.55,.7),textcoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#c8262d'},color='#a9242b',fontsize=11)
fig.suptitle('立型MC：安定限界とびびり振動形状\n基礎ばね支持・XYクロスFRF・4刃フルスロット零次近似',fontsize=13);fig.tight_layout();buf=io.BytesIO();fig.savefig(buf,format='png');plt.close(fig);buf.seek(0);left=Image.open(buf).convert('RGB').resize((1080,1440));left.save(R/'milling_selected.png')
m=build('foundation');grid,ids=vt.make_grid(type('M',(),{'parts':m.parts,'springs':[]})());grid.point_data['node_id']=ids;grid.cell_data['component']=m.element_parts;nodeidx=np.array(ids)-1
font=lambda n:ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),n)
colors=['#b5bbc3','#4f94c3','#a8cddc','#c93338','#e3ad46']

def make_scene(field, title, subtitle, detail, arrows=None):
 g=grid.copy();g.point_data['array_index']=nodeidx;p=pv.Plotter(off_screen=True,window_size=(920,480));p.set_background('white');actors=[];scale=65/max(np.linalg.norm(abs(field[nodeidx]),axis=1).max(),1e-30)
 for i,color in enumerate(colors):
  surf=g.extract_cells(np.where(np.array(m.element_parts)==i)[0]).extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface');ind=surf.point_data['array_index'].astype(int);points=surf.points.copy();du=field[ind]
  p.add_mesh(surf.copy(),style='wireframe',color='gray',opacity=.1);moving=surf.copy();p.add_mesh(moving,color=color);actors.append((moving,points,du))
 p.add_mesh(pv.Plane(center=(500,450,-70),direction=(0,0,1),i_size=1150,j_size=1050),color='#d9e4cd',opacity=.6)
 marks=[]
 if arrows is not None:
  for sign,color,index in [(1,'#c7262d',int(basis['work_index'])),(-1,'#247c3c',int(basis['tool_index']))]:
   point=pv.PolyData(CONTACT[None,:]);p.add_mesh(point,color=color,point_size=11,render_points_as_spheres=True);line=pv.Line(CONTACT,CONTACT+[50,0,0]);p.add_mesh(line,color=color,line_width=4);marks.append((point,line,sign,index))
 p.camera_position=[(2100,-1800,1600),(500,450,450),(0,0,1)];p.enable_parallel_projection();p.reset_camera();p.camera.zoom(.82);p.add_axes()
 def frame(phase):
  phasefactor=np.exp(1j*phase)
  for mesh,points,du in actors:mesh.points=points+scale*np.real(du*phasefactor)
  for point,line,sign,index in marks:
   pos=CONTACT+scale*np.real(field[index]*phasefactor);force=np.r_[np.real(arrows*phasefactor),0.]*sign;point.points=pos[None,:];line.points=np.array([pos,pos+100*force])
  p.render();im=Image.fromarray(p.screenshot(return_img=True));dr=ImageDraw.Draw(im);dr.rectangle((0,0,920,77),fill='white');dr.text((12,4),title,font=font(20),fill='black');dr.text((12,31),subtitle,font=font(18),fill='black');dr.text((12,56),f'{detail}　位相 {phase*180/np.pi:.0f}°',font=font(15),fill='black');dr.rectangle((0,450,920,480),fill='white');dr.text((12,455),'変位は各図で独立正規化。形状の振幅・変形倍率を図同士で比較しない。',font=font(15),fill='black');return im
 return p,frame

basis_frames=[]
for item in top:
 j=item['mode_number']-1;phi=basis['modes'][:,:,j];direction='XYZ'[np.argmax(abs(phi[int(basis['tool_index'])]))]
 p,frame=make_scene(phi,f"構造固有モード {j+1}：fn={item['natural_frequency_hz']:.1f} Hz",f"工具{direction}方向曲げが大きい構造モード／相対変位寄与指標 {item['relative_amplitude_share']*100:.1f}%",'びびり形状への寄与。単独の固有モードを表示')
 im=frame(0);p.close();im.save(R/f'structural_mode_{j+1}.png');basis_frames.append(im)
force=np.array(report['force_xy_real'])+1j*np.array(report['force_xy_imag'])
p,frame=make_scene(U,f"境界のびびり形状：fch={c['chatter_hz']:.1f} Hz",f"n={c['rpm']:.1f} min⁻¹、ap={c['depth_mm']:.3f} mm（選択点をモード基底で再計算）",'複数構造モード＋切削再生ループの複素形状',force)
frames=[]
for phase in np.linspace(0,2*np.pi,40,endpoint=False):
 im=Image.new('RGB',(2000,1440),'white');im.paste(left,(0,0));im.paste(frame(phase),(1080,0));im.paste(basis_frames[0],(1080,480));im.paste(basis_frames[1],(1080,960));frames.append(im)
p.close();frames[0].save(R/'milling_with_modes.png');frames[0].save(R/'milling_with_modes.gif',save_all=True,append_images=frames[1:],duration=80,loop=0)
assert ImageChops.difference(frames[0].crop((1200,85,1950,440)),frames[20].crop((1200,85,1950,440))).getbbox()
# Selected-condition global fields for ParaView and a complete machine-readable report.
for speed,r,U in states:
 g=grid.copy();g.point_data['U_real']=U[nodeidx].real;g.point_data['U_imag']=U[nodeidx].imag;g.save(R/f'chatter_mode_{speed}.vtu')
checks={'modal_vs_solver_FRF_relative_errors':errors,'selected_conditions':[{**r['condition'],'loop_residual':r['loop_relative_residual'],'modal_reconstruction_error':r['relative_reconstruction_error'],'dominant_modes':r['dominant_modes']} for _,r,_ in states],'display':'Right top: feedback boundary shape. Right middle/bottom: individual FE eigenmodes.'}
(R/'chatter_modes_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
