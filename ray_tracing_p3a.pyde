# This is the provided code for the ray tracing project.

# Philip Tarrer 
# project 3B

# The most important part of this code is the command interpreter, which
# parses the scene description (.cli) files.

from __future__ import division
import traceback

debug_flag = False   # print debug information when this is True

def setup():
    size(320, 320) 
    noStroke()
    colorMode(RGB, 1.0)  # Processing color values will be in [0, 1]  (not 255)
    background(0, 0, 0)
    frameRate(30)

# make sure proper error messages get reported when handling key presses
def keyPressed():
    try:
        handleKeyPressed()
    except Exception:
        traceback.print_exc()

# read and interpret a scene description .cli file based on which key has been pressed
def handleKeyPressed():
    if key == '1':
        interpreter("01_one_sphere.cli")
    elif key == '2':
        interpreter("02_three_spheres.cli")
    elif key == '3':
        interpreter("03_shiny_sphere.cli")
    elif key == '4':
        interpreter("04_many_spheres.cli")
    elif key == '5':
        interpreter("05_one_triangle.cli")
    elif key == '6':
        interpreter("06_icosahedron_and_sphere.cli")
    elif key == '7':
        interpreter("07_colorful_lights.cli")
    elif key == '8':
        interpreter("08_reflective_sphere.cli")
    elif key == '9':
        interpreter("09_mirror_spheres.cli")
    elif key == '0':
        interpreter("10_reflections_in_reflections.cli")
    elif key == '-':
        interpreter("11_star.cli")

lights = []
shapes = []
vertices = []
u = []
v = []
w = []


backgroundc = (0, 0, 0)
eyepos = (0, 0, 0)
fov = 0
material = None




# You should add code for each command that calls routines that you write.
# Some of the commands will not be used until Part B of this project.
def interpreter(fname):
    global lights, shapes, vertices, u, v, w, backgroundc, eyepos, fov, material
    
    reset_scene()  # you should initialize any data structures that you will use here
    
    fname = "data/" + fname
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()

    # parse the lines in the file in turn
    for line in lines:
        words = line.split()  # split up the line into individual tokens
        if len(words) == 0:   # skip empty lines
            continue
        if words[0] == 'sphere':
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            radius = float(words[1])
            shapes.append(Sphere(x, y, z, radius, material))
            
        elif words[0] == 'fov':
            fov = float(words[1])

        elif words[0] == 'eye':
            eyepos = [float(words[1]), float(words[2]), float(words[3])]
            
        elif words[0] == 'uvw':
            u.append([float(words[1]), float(words[2]), float(words[3])])
            v.append([float(words[4]), float(words[5]), float(words[6])])
            w.append([float(words[7]), float(words[8]), float(words[9])])

        elif words[0] == 'background':
            backgroundc = (float(words[1]), float(words[2]), float(words[3]))
            
        elif words[0] == 'light':
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            pos = [x, y, z]
            r = float(words[4])
            g = float(words[5])
            b = float(words[6])
            color = [r, g, b]
            lights.append(Light(pos, color))

        elif words[0] == 'surface':
            diffuse = [float(words[1]), float(words[2]), float(words[3])]
            ambient = [float(words[4]), float(words[5]), float(words[6])]
            specular = [float(words[7]), float(words[8]), float(words[9])]
            specpower = float(words[10])
            krefl = float(words[11])
            material = Material(diffuse, ambient, specular, specpower, krefl)

        elif words[0] == 'begin':
            vertices = []

        elif words[0] == 'vertex':
            vertices.append([float(words[1]), float(words[2]), float(words[3])])
            print "vertex: ", vertices[-1]

        elif words[0] == 'end':
            ab = [vertices[1][0] - vertices[0][0], vertices[1][1] - vertices[0][1], vertices[1][2] - vertices[0][2]]
            bc = [vertices[2][0] - vertices[1][0], vertices[2][1] - vertices[1][1], vertices[2][2] - vertices[1][2]]
            n = cross(ab, bc)
            shapes.append(Triangle(vertices[0], vertices[1], vertices[2], material, n))

        elif words[0] == 'render':
            render_scene()    # render the scene (this is where most of the work happens)
        elif words[0] == '#':
            pass  # ignore lines that start with the comment symbol (pound-sign)
        else:
            print ("unknown command: " + word[0])

# render the ray tracing scene
def render_scene():
    global lights, shapes, vertices, u, v, w, backgroundc, eyepos, fov, material

    for j in range(height):
        for i in range(width):
            d = 1 / tan(radians(fov / 2))
            U = ((2 * i) / width) - 1
            V = -1 * (((2 * j) / height) - 1)
            x = [-d * w[0][0], -d * w[0][1], -d * w[0][2]]
            y = [U * u[0][0], U * u[0][1], U * u[0][2]]
            z = [V * v[0][0], V * v[0][1], V * v[0][2]]
            direction = [x[0] + y[0] + z[0], x[1] + y[1] + z[1], x[2] + y[2] + z[2]]
            direction = normalize(direction)
            ray = Ray(eyepos, direction)
            hit = intersection(ray)
            if hit.shape is not None:
                c = color_scene(hit, ray, 10)
                c = color(c[0], c[1], c[2])
                set(i, j, c)
            else:
                c = color(backgroundc[0], backgroundc[1], backgroundc[2])
                set(i, j, c)
                       
def color_scene(hit, ray, depth):
    global lights, shapes, vertices, u, v, w, backgroundc, eyepos, fov, material
    tr = 0.0
    tg = 0.0
    tb = 0.0
    if (depth > 0 and hit.shape.material.kref > 0):
        N = hit.normal
        D = ray.direction
        DN = dot(N, scale(D, -1))
        scalar = 2 * DN
        temp = scale(N, scalar)
        R = add(temp, D)
        origin = add(hit.intersection, scale(N, 0.0001))
        newRay = Ray(origin, R)
        newHit = intersection(newRay)
        if newHit.shape is not None:
            newColor = color_scene(newHit, newRay, depth - 1)
            tr += newColor[0] * hit.shape.material.kref
            tg += newColor[1] * hit.shape.material.kref
            tb += newColor[2] * hit.shape.material.kref
        else:
            tr += backgroundc[0] * hit.shape.material.kref
            tg += backgroundc[1] * hit.shape.material.kref
            tb += backgroundc[2] * hit.shape.material.kref

    for light in lights:
        s = shadow(hit, light)
        L = sub(light.position, hit.intersection)
        L = normalize(L)
        N = hit.normal 
        D = hit.ray.direction
        temp = sub(L, D)
        H = normalize(temp)
        P = hit.shape.material.specpower

        specco = max(0, dot(H, N)) ** P
        diffuseco = max(0, dot(N, L))

        dr = hit.shape.material.diffuse[0] * light.color[0] * diffuseco
        dg = hit.shape.material.diffuse[1] * light.color[1] * diffuseco
        db = hit.shape.material.diffuse[2] * light.color[2] * diffuseco

        sr = hit.shape.material.specular[0] * light.color[0] * specco
        sg = hit.shape.material.specular[1] * light.color[1] * specco
        sb = hit.shape.material.specular[2] * light.color[2] * specco

        speccolor = [sr, sg, sb]
        diffusecolor = [dr, dg, db]
        temp = add(speccolor, diffusecolor)
        tr += temp[0] * s
        tg += temp[1] * s
        tb += temp[2] * s
    tr += hit.shape.material.ambient[0]
    tg += hit.shape.material.ambient[1]
    tb += hit.shape.material.ambient[2]
    return tr, tg, tb

# here you should reset any data structures that you will use for your scene (e.g. list of spheres)
def reset_scene():
    global lights, shapes, vertices, u, v, w, backgroundc, eyepos, fov, material
    lights = []
    shapes = []
    vertices = []
    u = []
    v = []
    w = []
    backgroundc = (0, 0, 0)
    eyepos = (0, 0, 0)
    fov = 0
    material = None

# prints mouse location clicks, for help debugging
def mousePressed():
    print ("You pressed the mouse at " + str(mouseX) + " " + str(mouseY))

# this function should remain empty for this assignment
def draw():
    pass

def intersection(ray):
    global shapes
    minT = 100000
    closestShape = None
    point = [0, 0, 0]
    normal = [0, 0, 0]
    for shape in shapes: 
        temp = shape.intersect(ray)
        if type(shape) is Sphere:
            if temp < minT:
                minT = temp
                closestShape = shape
                xt = ray.origin[0] + (ray.direction[0] * minT)
                yt = ray.origin[1] + (ray.direction[1] * minT)
                zt = ray.origin[2] + (ray.direction[2] * minT)
                point = [xt, yt, zt]
                temp = sub(point, [shape.x, shape.y, shape.z])
                normal = normalize(temp)
        else :
            if temp < minT:
                px = ray.origin[0] + (ray.direction[0] * temp)
                py = ray.origin[1] + (ray.direction[1] * temp)
                pz = ray.origin[2] + (ray.direction[2] * temp)
                point = [px, py, pz]
                ap = sub(point, shape.v1)
                ab = sub(shape.v2, shape.v1)
                if dot(shape.surfacenormal, ray.direction) < 0:
                    norm = shape.surfacenormal
                else :
                    norm = scale(shape.surfacenormal, -1)

                triple1 = dot(cross(ap, ab), norm)
                bp = sub(point, shape.v2)
                bc = sub(shape.v3, shape.v2)
                triple2 = dot(cross(bp, bc), norm)
                cp = sub(point, shape.v3)
                ca = sub(shape.v1, shape.v3)
                triple3 = dot(cross(cp, ca), norm)
                if (triple1 > 0) == (triple2 > 0) == (triple3 > 0):
                    if temp < minT:
                        minT = temp
                        closestShape = shape
                        point = [px, py, pz]
                        normal = norm
    return Hit(closestShape, normal, minT, point, ray)
       

def shadow(hit, light) :
    global shapes
    offset = scale(hit.normal, 0.0001)
    shadoworigin = add(hit.intersection, offset)
    shadowdirection = normalize(sub(light.position, hit.intersection))
    shadowray = Ray(shadoworigin, shadowdirection)
    
    temp = intersection(shadowray)
    if temp.shape is not None:
        if temp.t < distance(light.position, hit.intersection):
            return 0
    return 1

            


class Material(object) :
    def __init__(self, diffuse, ambient, specular, specpower, kref):
        self.diffuse = diffuse
        self.specular = specular
        self.ambient = ambient
        self.specpower = specpower
        self.kref = kref

class Sphere(object) :
    def __init__(self, x, y, z, radius, material):
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius
        self.material = material
    
    def intersect(self, ray):
        ox = ray.origin[0]
        oy = ray.origin[1]
        oz = ray.origin[2]
        dx = ray.direction[0]
        dy = ray.direction[1]
        dz = ray.direction[2]
        ux = ox - self.x
        uy = oy - self.y
        uz = oz - self.z

        a = dx ** 2 + dy ** 2 + dz ** 2
        b = 2 * (dx * ux + dy * uy + dz * uz)
        c = ux ** 2 + uy ** 2 + uz ** 2 - self.radius ** 2
        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return 100000
        else:
            t1 = (-b - sqrt(discriminant)) / (2 * a)
            t2 = (-b + sqrt(discriminant)) / (2 * a)
            if t1 < 0 and t2 < 0:
                return 100000
            elif t1 < 0:
                return t2
            elif t2 < 0:
                return t1
            else:
                return min(t1, t2)
            

class Light(object) :
    def __init__(self, position, color):
        self.position = position
        self.color = color

class Ray(object) :
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

class Hit(object) :
    def __init__(self, shape, normal, t, intersection, ray):
        self.shape = shape
        self.normal = normal
        self.t = t
        self.intersection = intersection
        self.ray = ray

class Triangle(object) :
    def __init__(self, v1, v2, v3, material, surfacenormal):
        self.surfacenormal = normalize(surfacenormal)
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.material = material

    def intersect(self, ray):
        denom = dot(self.surfacenormal, ray.direction)
        if denom != 0:
            t = dot(self.surfacenormal, (sub(self.v1, ray.origin))) / denom
            if t > 0:
                return t
            else:
                return 100000
        else : 
            return 100000

def normalize(vector):
    magnitude = sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    return [vector[0] / magnitude, vector[1] / magnitude, vector[2] / magnitude]
    

def mag(vector):
    return sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)

def dot(vector1, vector2):
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]

def cross(vector1, vector2):
    return [vector1[1] * vector2[2] - vector1[2] * vector2[1],
            vector1[2] * vector2[0] - vector1[0] * vector2[2],
            vector1[0] * vector2[1] - vector1[1] * vector2[0]]

def sub(vector1, vector2):
    return [vector1[0] - vector2[0], vector1[1] - vector2[1], vector1[2] - vector2[2]]

def add(vector1, vector2):
    return [vector1[0] + vector2[0], vector1[1] + vector2[1], vector1[2] + vector2[2]]

def scale(vector, scalar):
    return [vector[0] * scalar, vector[1] * scalar, vector[2] * scalar]

def distance(vector1, vector2):
    return sqrt((vector1[0] - vector2[0]) ** 2 + (vector1[1] - vector2[1]) ** 2 + (vector1[2] - vector2[2]) ** 2)
