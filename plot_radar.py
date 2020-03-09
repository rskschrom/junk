import numpy as np
from vispy import app, gloo
from get_radar import get_sweep_data

# radar colormaps
fil_zh = open('/home/robert/research/nws/zh2_map.rgb')
cdata = np.genfromtxt(fil_zh, skip_header=2)
zh_map = cdata/255

# read radar data
rpath = '/home/robert/research/nws/hires_loop/KICT20120414_223311_V06'
sweep_data = get_sweep_data(rpath)
zh = sweep_data[0]
ran = sweep_data[4]/1.e5
azi = sweep_data[5]*np.pi/180.

# set up radar grid
num_ran = len(ran)
num_azi = len(azi)
azmin = np.min(azi)
azmax = np.max(azi)
rmin = np.min(ran)
rmax = np.max(ran)

gspace = (rmax-rmin)/(num_ran-1)
azispace = 2.*np.pi/(num_azi-1)
azitmp = np.concatenate(([azi[-1]], azi, [azi[0]]))
azi_left = 0.5*(azitmp[:-2]+azitmp[1:-1])
azi_right = 0.5*(azitmp[2:]+azitmp[1:-1])

azispace_var = np.abs(azi_right-azi_left)
#azispace_var = 0.5*(np.abs(azitmp[1:-1]-azitmp[2:])+np.abs(azitmp[1:-1]-azitmp[:-2]))
#azispace_var = 0.5*(np.abs(azitmp[2:]-azitmp[:-2]))
print(azispace_var.shape, azi.shape)

#azispace_var = np.concatenate(([azispace_var[-1]], azispace_var, [azispace_var[0]]))
azispace_var[azispace_var>np.pi] = azispace_var[azispace_var>np.pi]-2.*np.pi

# quantize azimuth grid
dsazi = np.sign(azi[1]-azi[0])
#azi = azi[0]+dsazi*np.arange(num_azi)*azispace

th2, r2 = np.meshgrid(azi, ran, indexing='ij')
azispace2,_ = np.meshgrid(azispace_var, ran, indexing='ij')
thf = th2.flatten()
asf = azispace2.flatten()
rf = r2.flatten()

# set numpy structure arrays
nvert = 6
data = np.zeros(num_ran*num_azi*nvert, dtype=[('position', np.float32, 2),
                                              ('color', np.float32, 3)])
# polar grid (vertices of two triangles)
r1 = rf-gspace/2.
r2 = rf+gspace/2.
thet1 = thf-asf/2.
thet2 = thf+asf/2.
xpoints = np.stack((r1*np.cos(thet2),r2*np.cos(thet2),
                    r1*np.cos(thet1),r1*np.cos(thet1),
                    r2*np.cos(thet1),r2*np.cos(thet2)), axis=1).flatten()
ypoints = np.stack((r1*np.sin(thet2),r2*np.sin(thet2),
                    r1*np.sin(thet1),r1*np.sin(thet1),
                    r2*np.sin(thet1),r2*np.sin(thet2)), axis=1).flatten()
data['position'] = np.stack((xpoints, ypoints), axis=1)

# normalize field
field = zh[:]
fmin = 0.
fmax = 80.
field = (field-fmin)/(fmax-fmin)
field[field.mask] = 0.
field = np.clip(field, 0., 1.)

field_flat = field.flatten()

# set colors of field (actually plot data)
rc = 1.-np.zeros([num_ran*num_azi,1])
gc = 1.-np.zeros([num_ran*num_azi,1])
bc = 1.-np.zeros([num_ran*num_azi,1])

cidx = (field_flat[field_flat>0.]*255).astype(int)
rc[field_flat>0.,0] = zh_map[cidx,0]
gc[field_flat>0.,0] = zh_map[cidx,1]
bc[field_flat>0.,0] = zh_map[cidx,2]

rc = np.tile(rc, (1,nvert)).flatten()
gc = np.tile(gc, (1,nvert)).flatten()
bc = np.tile(bc, (1,nvert)).flatten()
data['color'] = np.stack((rc, gc, bc), axis=1)

c = app.Canvas(keys='interactive', size=(800,800))

vertex = """
uniform float scale;
uniform vec2 rel_pos;
attribute vec3 color;
attribute vec2 position;
varying vec3 v_color;
void main()
{
    gl_Position = vec4((position+rel_pos)*scale, 0.0, 1.0);
    v_color = color;
}
"""

fragment = """
varying vec3 v_color;
void main()
{
    gl_FragColor.rgb = v_color;
}
"""

program = gloo.Program(vertex, fragment, count=num_ran*num_azi*nvert)
program['color'] = data['color']
program['position'] = data['position']
program['rel_pos']    = (0.0, 0.0)
program['scale'] = 1.0

@c.connect
def on_resize(event):
    gloo.set_viewport(0, 0, *event.size)

@c.connect
def on_draw(event):
    gloo.clear((1,1,1,1))
    program.draw('triangles')

@c.connect
def on_mouse_move(event):
    win_size_x = c.size[0]
    win_size_y = c.size[1]

    rel_pos_init = program['rel_pos']
    sc = program['scale']

    if event.is_dragging and event.buttons[0] == 1:
        x0, y0 = event.last_event.pos[0], event.last_event.pos[1]
        x1, y1 = event.pos[0], event.pos[1]    

        norml_x0 = (2.0*x0/float(win_size_x))-1
        norml_y0 = (2.0*y0/float(win_size_y))-1
        norml_x1 = (2.0*x1/float(win_size_x))-1
        norml_y1 = (2.0*y1/float(win_size_y))-1

        program['rel_pos'] = (rel_pos_init[0]+(norml_x1-norml_x0)/sc, 
                              rel_pos_init[1]-(norml_y1-norml_y0)/sc)
        c.update()

@c.connect
def on_mouse_wheel(event):
    win_size_x = c.size[0]
    win_size_y = c.size[1]

    cur_pos_x = event.pos[0]
    cur_pos_y = event.pos[1]

    norm_x_pos = (2.0*cur_pos_x/win_size_x)-1.0
    norm_y_pos = (2.0*cur_pos_y/win_size_y)-1.0

    sign_scroll = event.delta[1]/abs(event.delta[1])

    program['scale'] = program['scale']*1.1**sign_scroll
    sc = program['scale'][0]

    program['rel_pos'] = program['rel_pos']+(-norm_x_pos*0.1*sign_scroll/sc, norm_y_pos*0.1*sign_scroll/sc)
    c.update()

c.show()
app.run();

