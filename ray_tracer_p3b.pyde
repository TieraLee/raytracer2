#X Tlee320
# This is the starter code for the CS 3451 Ray Tracing project.
#
# The most important part of this code is the interpreter, which will
# help you parse the scene description (.cli) files.
import math
from objects import *

scene = None
def setup():
    size(500, 500) 
    noStroke()
    colorMode(RGB, 1.0)  # Processing color values will be in [0, 1]  (not 255)
    background(0, 0, 0)

# read and interpret the appropriate scene description .cli file based on key press
def keyPressed():  
    if key == '1':
        interpreter("i1.cli")
    elif key == '2':
        interpreter("i2.cli")
    elif key == '3':
        interpreter("i3.cli")
    elif key == '4':
        interpreter("i4.cli")
    elif key == '5':
        interpreter("i5.cli")
    elif key == '6':
        interpreter("i6.cli")
    elif key == '7':
        interpreter("i7.cli")
    elif key == '8':
        interpreter("i8.cli")
    elif key == '9':
        interpreter("i9.cli")
    elif key == '0':
        interpreter("i10.cli")


def interpreter(fname):
    global scene
    scene = Scene()  
    fname = "data/" + fname
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()

    # parse each line in the file in turn
    for line in lines:
        words = line.split()  # split the line into individual tokens
        if len(words) == 0:   # skip empty lines
            continue
        if words[0] == 'sphere':
            radius = float(words[1])
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            #Add sphere to object list
            newSphere = sphereObject(x,y,z,radius,scene.current_surface)
            # m_sphere.set_surface(scene.current_surface)
            scene.sceneObjects.append(newSphere)
            #x,y,z are location , r is how big
            # call your sphere creation routine here
            # for example: create_sphere(radius,x,y,z)
        elif words[0] == 'fov':
            global fov 
            fov = float(words[1])
            scene.fov = fov
        elif words[0] == 'background':
            r = float(words[1])
            g = float(words[2])
            b = float(words[3])
            scene.backgroundColor = (r,g,b)
        elif words[0] == 'light':
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            r = float(words[4])
            g = float(words[5])
            b = float(words[6])
            scene.lightSources.append(Light(x,y,z,r,g,b))
            
        elif words[0] == 'surface':
            dif_r = float(words[1])
            dif_g = float(words[2])
            dif_b = float(words[3])
            amb_r = float(words[4])
            amb_g = float(words[5])
            amb_b = float(words[6])
            spec_r = float(words[7])
            spec_g = float(words[8])
            spec_b = float(words[9])
            spec_power = float(words[10])
            reflc = float(words[11])
            scene.current_surface = dif_r, dif_g, dif_b,amb_r, amb_g, amb_b, spec_r, spec_g, spec_b, spec_power, reflc
            
        elif words[0] == 'begin':
            scene.begin = True

        elif words[0] == 'vertex':
            if scene.begin:
                x = float(words[1])
                y = float(words[2])
                z = float(words[3])
                scene.vertices.append(Vector(x,y,z))
        elif words[0] == 'end':
            scene.end_vertex()

        elif words[0] == 'write':
            render_scene()    # render the scene
            save(words[1])  # write the image to a file
            pass

# render the ray tracing scene
def render_scene():
    global scene
    for j in range(height):
        for i in range(width):
            ray = Ray(0,0,0,i,j)
            ray.convert_to_per(scene.fov)
            hit = scene.calculate_intersection(ray,1)
            if hit.closest_object == None:
                r,g,b = scene.backgroundColor
                bg = color(r,g,b)
                set(i,j, bg)
            else:
                r, g ,b = hit.f_color
                bg = color(r, g, b)
                set(i,j, bg)
    

# should remain empty for this assignment
def draw():
    pass