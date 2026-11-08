from pathlib import Path
import sys,numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
from foundation_model import build
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent/'tools'))
import vtu_tools as vt
m=build();R=D/'foundation/results';g,ids=vt.make_grid(type('M',(),{'parts':m.parts,'springs':[]})());g.point_data['node_id']=ids
coords={c[0]:np.array(c[1:]) for part in m.parts for c in part.coords}
p=pv.Plotter(off_screen=True,window_size=(1500,1000));p.set_background('white');p.add_mesh(g,color='#b3bfcc',opacity=.12,show_edges=False)
p.add_mesh(pv.Plane(center=(500,450,-180),direction=(0,0,1),i_size=1150,j_size=1050),color='#a9bc98',opacity=.5)
network={'solids':g};points=[];labs=[]
for foot in m.foundations:
 center=np.array(foot['xyz']);ref=center+[0,0,-40];ground=center+[0,0,-180]
 faces=np.array([coords[n] for n in foot['surface_nodes']]);segs=np.array([x for v in faces for x in (ref,v)])
 fan=pv.line_segments_from_points(segs);p.add_mesh(fan,color='#8d2bb5',line_width=2);p.add_points(faces,color='#8d2bb5',point_size=6,render_points_as_spheres=True)
 spring=pv.Line(ref,ground);p.add_mesh(spring,color='#cb2834',line_width=6)
 p.add_points(np.array([ref]),color='#e3b420',point_size=14,render_points_as_spheres=True);p.add_points(np.array([ground]),color='#28753b',point_size=15,render_points_as_spheres=True)
 points.append(ground+[0,0,-25]);labs.append(f"{foot['name']}: B{foot['ref']} - G{foot['ground']}\n{len(faces)} surface nodes / {len(faces)*3} MPC eqs")
 network[foot['name']+'_MPC_display']=fan;network[foot['name']+'_spring_display']=spring
p.add_point_labels(points,labs,font_size=16,text_color='black',shape_color='white',shape_opacity=.9,always_visible=True,point_size=0)
p.camera_position=[(1800,-1600,-1250),(500,450,200),(0,0,1)];p.reset_camera();p.camera.zoom(.85);p.add_axes()
im=Image.fromarray(p.screenshot(return_img=True));p.close();dr=ImageDraw.Draw(im);font=lambda n:ImageFont.truetype(str(Path.home()/'.fonts/meiryo.ttc'),n)
dr.rectangle((0,0,1500,115),fill='white');dr.text((15,8),'追加モデル：ベッド底面の4取付面 → MPC → ベッド側代表点 → XYZばね → 固定地面側代表点',font=font(23),fill='black')
dr.text((15,48),'下側から見た図。紫：面節点とMPC　黄：ベッド側代表点　赤：ばね　緑：固定地面側代表点',font=font(21),fill='black')
dr.text((15,80),'各組の代表点の実距離は0 mm。表示だけBをz=-40、Gをz=-180 mmへ移動（赤線の表示長140 mm）。',font=font(20),fill='black')
dr.rectangle((0,920,1500,1000),fill='white');dr.text((15,930),'1か所あたり kX=kY=10,000 N/mm、kZ=100,000 N/mm（説明用の仮定値）。地面側だけXYZ固定。',font=font(21),fill='black')
dr.text((15,965),'ベッド底面全体の固定は解除。4か所以外の底面節点は直接固定しない。地盤のソリッド要素は作っていない。',font=font(20),fill='black')
im.save(R/'foundation_mpc.png');pv.MultiBlock(network).save(R/'foundation_mpc_display.vtm')
print([(f['name'],f['ref'],f['ground'],len(f['surface_nodes'])) for f in m.foundations])
