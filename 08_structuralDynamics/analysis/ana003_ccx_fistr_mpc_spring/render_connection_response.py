"""CAR1 actual-node visualization and forced-response GIF above its FRF peak."""
from pathlib import Path
import sys,io,json
import numpy as np
from PIL import Image,ImageDraw,ImageFont,ImageChops
import pyvista as pv
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'tools'))
import mpc_springs as ms
import vtu_tools as vt
from plot_modes import plt
from post_frequency import read_front
R=D/'results/frf'; R.mkdir(exist_ok=True)
font=lambda s:ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),s)
m=ms.build('rbe2');g,ids=vt.make_grid(type('M',(),{'parts':m.parts,'springs':[]})());g.point_data['node_id']=ids
surf=g.extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface')
coords={n:np.array([x,y,z]) for p in m.parts for n,x,y,z in p.coords};coords.update({n:np.array([x,y,z]) for n,x,y,z in m.extra_nodes})
a,b=22259,22261
assert np.array_equal(coords[a],coords[b]);assert any(s[1]==2 and s[3]==a and s[4]==b for s in m.springs)
node=pv.PolyData(np.array([coords[a],coords[b]]));node.point_data['node_id']=[a,b]
# Actual solver geometry, preserving coincident but distinct node IDs.
pv.MultiBlock({'solids':g,'CAR1_reference_nodes':node}).save(R/'car1_actual.vtm')
local=[]
for part in (0,1):
 s=g.extract_cells(np.where(g.cell_data['part']==part)[0])
 s=s.clip_box((245,375,105,195,185,285),invert=False,crinkle=True).extract_surface(nonlinear_subdivision=0,algorithm='dataset_surface')
 local.append(s)
lift=np.array([0.,0.,150.]); bed=local[0];table=local[1].copy();table.points+=lift
pa=coords[a];pb=coords[b]+lift
line=pv.Line(pa,pb);line.cell_data['element_id']=np.full(line.n_cells,4289)
refs=pv.PolyData(np.array([pa,pb]));refs.point_data['node_id']=[a,b]
pv.MultiBlock({'bed_rail':bed,'table_carriage_exploded_Z_plus150mm':table,'spring_4289_SP0':line,'reference_nodes':refs}).save(R/'car1_exploded.vtm')
def setup(p,animated=False):
 p.set_background('white');p.add_mesh(bed,color='#abb5c0',show_edges=True,edge_color='#647080',opacity=.8 if animated else .12)
 dra=pa if animated else pa+np.array([0,0,30.]);drb=pb if animated else pb-np.array([0,0,30.])
 p.add_mesh(table,style='wireframe',color='#4279ae',opacity=.18)
 moving=table.copy();p.add_mesh(moving,color='#63b3e5',show_edges=True,edge_color='#31577a',opacity=1. if animated else .12)
 spring=pv.Line(dra,drb);p.add_mesh(spring,color='#d52d34',line_width=7)
 pointmesh=pv.PolyData(np.array([dra,drb]));p.add_mesh(pointmesh,color='#eab600',point_size=16,render_points_as_spheres=True)
 p.add_point_labels(np.array([dra+[-35,0,0],drb+[15,0,0]]),['22259 : BED','22261 : TABLE'],font_size=17,point_size=0,always_visible=True,shape_color='white',shape_opacity=.9,text_color='black')
 for part,ref,offset in [(m.parts[0],dra,np.zeros(3)),(m.parts[1],drb,lift)]:
  points=ms._patch_nodes(part,(270,350,130,170),230)
  seg=[]
  for n,x,y in points:seg.extend([ref,coords[n]+offset])
  p.add_mesh(pv.line_segments_from_points(np.array(seg)),color='#8d2bb5',line_width=2,opacity=1.)
  if not animated:
   facepts=np.array([coords[n]+offset for n,x,y in points]);p.add_points(facepts,color='#8d2bb5',point_size=7,render_points_as_spheres=True)
   p.add_point_labels([ref+np.array([55,0,0])],[f'{len(points)} nodes / 111 equations'],font_size=14,point_size=0,always_visible=True,text_color='#6c188e',shape_color='white')
 p.camera_position=[(680,-380,550),(310,150,300),(0,0,1)];p.reset_camera();p.camera.zoom(.70)
 p.add_axes()
 return moving,spring,pointmesh
p=pv.Plotter(shape=(1,2),off_screen=True,window_size=(1400,780));p.subplot(0,0);p.set_background('white')
p.add_mesh(surf,color='#a6b8c9',opacity=.38,show_edges=False);p.add_points(pa,color='red',point_size=20,render_points_as_spheres=True)
p.add_point_labels([pa],['CAR1: 22259 + 22261'],font_size=17,always_visible=True,text_color='black',shape_color='white')
p.camera_position=[(2100,-1700,1400),(500,380,400),(0,0,1)];p.reset_camera();p.camera.zoom(.8)
p.subplot(0,1);setup(p)
def display_xy(point):
 p.renderer.SetWorldPoint(*point,1.);p.renderer.WorldToDisplay();x,y,z=p.renderer.GetDisplayPoint();return x,780-y
upper_xy=display_xy(pb-[0,0,30]);lower_xy=display_xy(pa+[0,0,30])
im=Image.fromarray(p.screenshot(return_img=True));p.close();dr=ImageDraw.Draw(im)
dr.rectangle((0,0,1400,90),fill='white');dr.text((15,8),'要素4289（SP0）が結ぶ点：22259（ベッド側）と22261（テーブル側）',font=font(24),fill='black')
dr.text((15,48),'左：実際の位置（2点は重なる）　　右：CAR1拡大・分解表示（テーブル側をZ方向へ150 mm移動）',font=font(20),fill='black')
dr.rectangle((0,680,1400,780),fill='white')
dr.text((15,690),'紫の点・線：各面37節点を代表点に結ぶMPC（1面111式）。黄色：代表点。赤：2節点ばね。',font=font(21),fill='black')
dr.text((15,728),'見やすく代表点も面から30 mm離して表示。実座標は両点とも (310,150,230) mm。線の長さ・方向は表示用。',font=font(19),fill='black')
# Actual coordinates and zero length must not be inferred from the exploded view.
dr.rounded_rectangle((16,100,670,200),radius=10,fill='white',outline='#c52222',width=2)
dr.text((30,110),'22259 ↔ 22261 の実距離：0 mm',font=font(27),fill='#c52222')
dr.text((30,151),'両節点の実座標：(310, 150, 230) mm',font=font(21),fill='black')
x1,y1=upper_xy;x2,y2=lower_xy;xm=min(x1,x2)-130
for x,y in [(x1,y1),(x2,y2)]: dr.line((x-10,y,xm,y),fill='#b92323',width=2)
dr.line((xm,y1,xm,y2),fill='#b92323',width=2)
for yy,sign in [(y1,1),(y2,-1)]:dr.polygon([(xm,yy),(xm-6,yy+sign*13),(xm+6,yy+sign*13)],fill='#b92323')
dr.rectangle((xm-75,(y1+y2)/2-30,xm+100,(y1+y2)/2+28),fill='white')
dr.text((xm-70,(y1+y2)/2-26),'表示上 90 mm',font=font(19),fill='#b92323')
dr.text((xm-70,(y1+y2)/2),'実距離 0 mm',font=font(19),fill='#b92323')
im.save(R/'car1_nodes.png')
network={'bed':bed,'table_exploded':table}
network['spring_4289_display']=pv.Line(pa+[0,0,30],pb-[0,0,30])
for part,ref,offset,name in [(m.parts[0],pa+[0,0,30],np.zeros(3),'BED'),(m.parts[1],pb-[0,0,30],lift,'TABLE')]:
 pts=ms._patch_nodes(part,(270,350,130,170),230);seg=[]
 for n,x,y in pts:seg.extend([ref,coords[n]+offset])
 network[name+'_MPC_display_links']=pv.line_segments_from_points(np.array(seg))
 nodes=pv.PolyData(np.array([coords[n]+offset for n,x,y in pts]));nodes.point_data['node_id']=[n for n,x,y in pts];network[name+'_surface_nodes']=nodes
pv.MultiBlock(network).save(R/'car1_mpc_network.vtm')
if '--diagram-only' in sys.argv: sys.exit(0)
# Reuse the whole-machine renderer for the response GIFs.
import runpy
runpy.run_path(str(D/'render_whole_peak.py'),run_name='__main__')
