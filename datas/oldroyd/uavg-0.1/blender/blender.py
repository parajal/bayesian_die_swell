import bpy, bmesh, math, numpy as np, sys

# ---------- 1. GEOMETRY ----------
x = np.loadtxt('curve4_x.txt')
x, r = np.loadtxt('lambda_1.5_beta_0.25')
n_theta = 64

verts, faces = [], []
th = np.linspace(0.0, 2.0*math.pi, n_theta, endpoint=False)
ct, st = np.cos(th), np.sin(th)
nx = len(x)
for i in range(nx):
    for j in range(n_theta):
        verts.append((float(x[i]), float(r[i]*ct[j]), float(r[i]*st[j])))
idx = lambda i, j: i*n_theta + (j % n_theta)
for i in range(nx-1):
    for j in range(n_theta):
        faces.append((idx(i,j), idx(i,j+1), idx(i+1,j+1), idx(i+1,j)))
# end caps
faces.append(tuple(idx(0, j) for j in range(n_theta-1, -1, -1)))
faces.append(tuple(idx(nx-1, j) for j in range(n_theta)))

# clean scene
bpy.ops.wm.read_factory_settings(use_empty=True)
me = bpy.data.meshes.new("swell")
me.from_pydata(verts, [], faces)
me.validate()
me.update()
for p in me.polygons:
    p.use_smooth = True
ob = bpy.data.objects.new("swell", me)
bpy.context.scene.collection.objects.link(ob)

# ---------- 3. MATERIAL ----------
mat = bpy.data.materials.new("extrudate")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.8, 0.2, 0.05, 1.0)
bsdf.inputs["Roughness"].default_value = 0.35
if "Metallic" in bsdf.inputs:
    bsdf.inputs["Metallic"].default_value = 0.0
ob.data.materials.append(mat)

# ---------- 2. CAMERA ----------
cam_data = bpy.data.cameras.new("cam")
cam_data.type = 'ORTHO'
cam_data.ortho_scale = 9.0
cam = bpy.data.objects.new("cam", cam_data)
cam.location = (1.0, -20.0, 0.0)
bpy.context.scene.collection.objects.link(cam)
# aim at (1.0, 0, 0): looking along +Y
cam.rotation_euler = (math.pi/2, 0.0, 0.0)
bpy.context.scene.camera = cam

# ---------- 4. LIGHT ----------
lt = bpy.data.lights.new("key", type='POINT')
lt.energy = 1000.0
lt.shadow_soft_size = 1.5
light = bpy.data.objects.new("key", lt)
light.location = (2.5, -10.0, 10.0)
bpy.context.scene.collection.objects.link(light)

# soft ambient so the silhouette reads
world = bpy.data.worlds.new("W")
bpy.context.scene.world = world
world.use_nodes = True
world.node_tree.nodes["Background"].inputs[0].default_value = (0.05, 0.055, 0.07, 1.0)
world.node_tree.nodes["Background"].inputs[1].default_value = 1.0

# ---------- 5. RENDER ----------
sc = bpy.context.scene
sc.render.resolution_x, sc.render.resolution_y = 1600, 1000
sc.render.resolution_percentage = 100
sc.render.filepath = "/tmp/blend/swell.png"
sc.render.image_settings.file_format = 'PNG'

engine = sys.argv[-1] if len(sys.argv) > 1 else 'CYCLES'
sc.render.engine = engine
if engine == 'CYCLES':
    sc.cycles.device = 'CPU'
    sc.cycles.samples = 48
    sc.cycles.use_denoising = True
    sc.cycles.max_bounces = 4
print("verts", len(verts), "faces", len(faces), "engine", sc.render.engine)
bpy.ops.render.render(write_still=True)
print("DONE")