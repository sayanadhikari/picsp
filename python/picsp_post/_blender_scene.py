"""Blender side of picsp_post.blender: build and render the exported scene.

Run by Blender, not by Python:  blender -b --factory-startup -P _blender_scene.py -- scene.json
"""
import json
import math
import sys

import bpy
from mathutils import Matrix, Vector

scene_file = sys.argv[sys.argv.index("--") + 1]
with open(scene_file) as f:
    S = json.load(f)

for ob in list(bpy.data.objects):
    bpy.data.objects.remove(ob, do_unlink=True)

sc = bpy.context.scene
sc.render.engine = "CYCLES"
sc.cycles.device = "CPU"
sc.cycles.samples = S["samples"]
sc.cycles.use_denoising = True
sc.cycles.max_bounces = 6
sc.render.film_transparent = True
sc.render.resolution_x, sc.render.resolution_y = S["width"], S["height"]
sc.render.resolution_percentage = 100
sc.render.image_settings.file_format = "PNG"
sc.render.image_settings.color_mode = "RGBA"
sc.render.image_settings.color_depth = "16"
sc.render.filepath = S["output"]
# data colours must come out as they went in: no filmic tone mapping
sc.view_settings.view_transform = "Standard"
sc.view_settings.look = "None"

# soft white studio light from the environment
world = bpy.data.worlds.new("studio")
sc.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
bg.inputs[0].default_value = (1, 1, 1, 1)
bg.inputs[1].default_value = 0.75


def material(name, kind, opacity):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    attr = nt.nodes.new("ShaderNodeVertexColor")
    attr.layer_name = "Col"
    if kind == "flat":
        em = nt.nodes.new("ShaderNodeEmission")
        em.inputs["Strength"].default_value = 1.0
        nt.links.new(attr.outputs["Color"], em.inputs["Color"])
        nt.links.new(em.outputs[0], out.inputs["Surface"])
        return m
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.inputs["Roughness"].default_value = 0.42
    bsdf.inputs["Specular IOR Level"].default_value = 0.35
    bsdf.inputs["Alpha"].default_value = opacity
    nt.links.new(attr.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs[0], out.inputs["Surface"])
    return m


for i, mesh in enumerate(S["meshes"]):
    bpy.ops.wm.ply_import(filepath=mesh["path"])
    ob = bpy.context.selected_objects[0]
    if ob.data.color_attributes:
        ob.data.color_attributes[0].name = "Col"
    ob.data.materials.append(material(f"m{i}", mesh["kind"], mesh["opacity"]))
    if mesh["lines"]:
        ob.data.shade_smooth()
    elif mesh["kind"] == "surface":
        # CAD and STL parts often have mixed triangle winding: orient the normals outwards
        bpy.context.view_layer.objects.active = ob
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.select_all(action="SELECT")
        if mesh["smooth"]:
            bpy.ops.mesh.remove_doubles(threshold=1e-6 * S["diag"])
        bpy.ops.mesh.normals_make_consistent(inside=False)
        bpy.ops.object.mode_set(mode="OBJECT")
        if mesh["smooth"]:
            # smooth curved surfaces but keep the sharp edges of CAD parts
            bpy.ops.object.shade_smooth_by_angle(angle=math.radians(35))
    if mesh["kind"] == "flat":
        # a data slice: no shadows on it and none cast by it
        ob.visible_shadow = False

# camera from the pyvista view
C = S["camera"]
pos, foc, up = Vector(C["position"]), Vector(C["focal_point"]), Vector(C["up"])
fwd = (foc - pos).normalized()
right = fwd.cross(up).normalized()
upv = right.cross(fwd).normalized()
rot = Matrix((right, upv, -fwd)).transposed()
cam_data = bpy.data.cameras.new("cam")
cam = bpy.data.objects.new("cam", cam_data)
sc.collection.objects.link(cam)
cam.matrix_world = Matrix.Translation(pos) @ rot.to_4x4()
cam_data.sensor_fit = "VERTICAL"
if C["parallel"]:
    cam_data.type = "ORTHO"
    cam_data.ortho_scale = 2 * C["parallel_scale"] * max(1.0, S["width"] / S["height"])
else:
    cam_data.angle = math.radians(C["view_angle"])
dist = (foc - pos).length
cam_data.clip_start = dist * 1e-3
cam_data.clip_end = dist * 20 + S["diag"] * 10
sc.camera = cam

# key light: a sun from the upper left behind the viewer, soft shadows
sun_data = bpy.data.lights.new("key", "SUN")
sun_data.energy = 2.2
sun_data.angle = math.radians(12)
sun = bpy.data.objects.new("key", sun_data)
sc.collection.objects.link(sun)
d = (-0.55 * right + 0.8 * upv - 0.6 * fwd).normalized()   # direction towards the light
z = d
x = z.cross(Vector((0, 0, 1)) if abs(z.z) < 0.99 else Vector((0, 1, 0))).normalized()
y = z.cross(x)
sun.matrix_world = Matrix((x, y, z)).transposed().to_4x4()
# a shadow cast by the studio light on a surface coloured by data (e.g. a photoemission map,
# which has its own physical shadow) would be read as data: no key-light shadows then
if any(m["kind"] == "surface" and m.get("data") for m in S["meshes"]):
    sun_data.use_shadow = False

bpy.ops.render.render(write_still=True)
