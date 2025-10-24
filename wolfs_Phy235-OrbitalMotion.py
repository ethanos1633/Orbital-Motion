from vpython import *
from math import sqrt, pi, cos, sin

# ------------------------------
# Scene setup
# ------------------------------
scene.title = "Earth-Mars Hohmann Transfer"
scene.background = color.white
scene.width = 800
scene.height = 600
scene.autoscale = 0
scene.range = 3e11
scene.center = vector(0,0,0)

# ------------------------------
# Constants
# ------------------------------
G = 6.7e-11  # gravitational constant

# ------------------------------
# Sun
# ------------------------------
sun = sphere(pos=vector(0,0,0), radius=10*0.7e9, color=color.yellow, mass=2e30)
sun.p = vector(0,0,0)

# ------------------------------
# Earth
# ------------------------------
earth_radius = 1.5e11
earth = sphere(pos=vector(earth_radius,0,0), radius=6.4e6*1000, color=color.blue, mass=6e24)
earth.v = vector(0, sqrt(G*sun.mass/mag(earth.pos)), 0)
earth.p = earth.mass * earth.v

# ------------------------------
# Mars
# ------------------------------
mars_radius = 2.28e11
r1 = earth_radius
r2 = mars_radius
t_trans = pi * sqrt(((r1 + r2)**3) / (8 * G * sun.mass))
n2 = sqrt(G * sun.mass / r2**3)
theta = pi - n2 * t_trans

mars_pos = vector(r2*cos(theta), r2*sin(theta), 0)
mars_v = vector(-sqrt(G*sun.mass/r2)*sin(theta), sqrt(G*sun.mass/r2)*cos(theta), 0)
mars = sphere(pos=mars_pos, radius=3.4e6*1000, color=color.red, mass=6.4e23)
mars.v = mars_v
mars.p = mars.mass * mars.v

# ------------------------------
# Rocket (start at Earth)
# ------------------------------
rocket = sphere(pos=earth.pos, radius=5e6, color=color.green, mass=1e4)
rocket.v = vector(0,0,0)
rocket.p = rocket.mass * rocket.v

# ------------------------------
# Orbit curves
# ------------------------------
for obj in [sun, earth, mars, rocket]:
    obj.orbit = curve(color=obj.color, radius=2e9)

# ------------------------------
# Graphs for x and y
# ------------------------------
def create_xy_graph(title, color_x, color_y):
    g = graph(width=500, height=250, title=title, xtitle="Time (s)", ytitle="Position (m)")
    gx = gcurve(graph=g, color=color_x, label="x")
    gy = gcurve(graph=g, color=color_y, label="y")
    return gx, gy

sun_x, sun_y = create_xy_graph("Sun x and y", color.blue, color.red)
earth_x, earth_y = create_xy_graph("Earth x and y", color.blue, color.red)
mars_x, mars_y = create_xy_graph("Mars x and y", color.orange, color.purple)
rocket_x, rocket_y = create_xy_graph("Rocket x and y", color.green, color.black)

# ------------------------------
# Simulation parameters
# ------------------------------
dt = 3600
time = 0
max_time = 5e7
launch_done = False

# ------------------------------
# Gravitational force
# ------------------------------
def grav_force(obj1, obj2):
    r_vec = obj1.pos - obj2.pos
    r_mag = mag(r_vec)
    if r_mag < 1e5:
        return vector(0,0,0)
    return -G * obj1.mass * obj2.mass * r_vec / r_mag**3

# ------------------------------
# Hohmann delta-v
# ------------------------------
def hohmann_delta_v(r1, r2):
    v1 = sqrt(G*sun.mass/r1)
    a = (r1 + r2)/2
    v_transfer = sqrt(G*sun.mass*(2/r1 - 1/a))
    return v_transfer - v1

delta_v = hohmann_delta_v(earth_radius, mars_radius)
print("Delta-v Earth->Mars: {:.2f} m/s".format(delta_v))

# ------------------------------
# Simulation loop
# ------------------------------
while time < max_time:
    rate(200)
    
    # Forces
    F_SE = grav_force(earth, sun)
    F_SM = grav_force(mars, sun)
    F_RS = grav_force(rocket, sun)

    # Update planet momenta
    earth.p += F_SE*dt
    mars.p += F_SM*dt

    # Launch rocket once
    if not launch_done:
        rocket.v = earth.v + norm(earth.v)*delta_v
        rocket.p = rocket.mass * rocket.v
        launch_done = True
        launch_time = time

    # Update rocket momentum
    rocket.p += F_RS*dt

    # Update positions
    for obj in [earth, mars, rocket]:
        obj.pos += obj.p/obj.mass * dt
        obj.orbit.append(pos=obj.pos)

    # Update graphs
    sun_x.plot(pos=(time, sun.pos.x))
    sun_y.plot(pos=(time, sun.pos.y))
    earth_x.plot(pos=(time, earth.pos.x))
    earth_y.plot(pos=(time, earth.pos.y))
    mars_x.plot(pos=(time, mars.pos.x))
    mars_y.plot(pos=(time, mars.pos.y))
    rocket_x.plot(pos=(time, rocket.pos.x))
    rocket_y.plot(pos=(time, rocket.pos.y))

    # Arrival check
    if mag(rocket.pos - mars.pos) < 1e9 and 'arrival_time' not in globals():
        arrival_time = time - launch_time
        print("Rocket reached Mars orbit after {:.2f} days.".format(arrival_time/86400))

    time += dt
