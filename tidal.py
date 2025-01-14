from __future__ import division
from visual import *
from visual.graph import *
#Constants------------------------------------------------------------

G = 6.67e-11
Water_Mass = 2.99e-23
Moon_Mass = 7.34e22
Earth_Mass = 5.97e24
Earth_Radius = 6371
GMm = G*Water_Mass*Moon_Mass
#---------------------------------------------------------------------
Scene = display(forward=vector(-.1,-.25,-1))
Moon = sphere(pos= (50000,0,0), radius = 1737)
Earth = sphere(pos = (0,0,0),radius =
Earth_Radius,material=materials.earth,opacity=.7)
Earth.i = 8.034e37
Earth.o = 7.27e-5
Earth.p = Earth.o*Earth.i
#---------------------------------------------------------------------
Ocean = [] #list of water molecules
Ocean_Frame = frame() #frame that rotates with Phi
#Create the Ocean

for phi in arange(0,pi,pi/20):
    for theta in arange(0,2*pi,pi/20):
        WaterMol = sphere(frame = Ocean_Frame,
pos=(6400*cos(theta)*sin(phi),6400*cos(phi),
6400*sin(theta)*sin(phi)),radius=100, color=color.cyan)
    WaterMol.r = sqrt(WaterMol.pos.x*WaterMol.pos.x +
    WaterMol.pos.z*WaterMol.pos.z)
    WaterMol.beta = theta
    Ocean.append(WaterMol)
#Time/Angle Parameters
#------------------------------------------------------------

Theta=0 #angle of moon orbit
dTheta=.001
t = 0 #time in seconds
dt = 29.5*24*60*60/(2000*3.14)
Phi = 0 #angle of earth rotation
dPhi = Earth.o*dt
#---------------------------------------------------------------------
#This function updates the position of the water molecule baded on its
#fixed radius and changing angle
def update_pos( sphere, alpha, radius ):
    sphere.pos = (radius*cos(alpha),sphere.pos.y,-radius*sin(alpha))
#---------------------------------------------------------------------

#This function determines the real position of an object in a frame

def world_space_pos(frame, local):
    x_axis = norm(frame.axis)
    z_axis = norm(cross(frame.axis, frame.up))
    y_axis = norm(cross(z_axis, x_axis))
    return frame.pos+local.x*x_axis+local.y*y_axis+local.z*z_axis
#Plot used for tests
#-----------------------------------------------------------------
display1 = gdisplay(x=400,y=0,width=700,height=450,xtitle="Angle",ytitle="Years")
K = gcurve(gdisplay=display1,color=color.blue)
#---------------------------------------------------------------------

while Theta < pi/120:
    Moon.rpos = vector(384400*cos(Theta),0,-384400*sin(Theta))=#updates real position of the moon
    Moon.opp = vector( -Moon.rpos.x, 0, -Moon.rpos.z) #opposite position of moon, used for gravitational field
    Moon.pos =(50000*cos(Theta),0,-50000*sin(Theta)) #updates simulated position of the moon
    Earth.rotate(angle=dPhi,axis=(0,1,0)) #rotates the earth by dPhi
    Ocean_Frame.rotate(angle=dPhi,axis=(0,1,0)) #rotates the ocean frame by dPhi

for object in Ocean: #determines Fg for each water molecule and updates its position accordingly
    object.rpos = world_space_pos(Ocean_Frame,object.pos)
    while object.beta >= 2*pi:
        object.beta = object.beta - 2*pi
    if object.beta > Theta and object.beta < (Theta + pi/2):
        object.Fg = -GMm/mag2(Moon.rpos-object.rpos)
    elif object.beta > (Theta + pi/2) and object.beta < (Theta +
    pi):
        object.Fg = GMm/mag2(Moon.opp-object.rpos)
    elif object.beta > (Theta + pi) and object.beta < (Theta +
    3*pi/2):
        object.Fg = -GMm/mag2(Moon.opp-object.rpos)
    else:
        object.Fg = GMm/mag2(Moon.rpos-object.rpos)
        object.Fg = object.Fg*abs(sin(object.beta)*cos(object.beta))
        object.beta = object.beta + object.Fg/Water_Mass/Earth_Radius
        + dPhi
        update_pos(object,object.beta,object.r)

    Theta=Theta+dTheta
    t=t+dt
    Phi = Phi + dPhi

    rate(10)
    #---------------------------------------------------------------------
#Computes the ratio of torques by high tide and low tide

High_Tide_Grav = 0
Low_Tide_Grav = 0
for object in Ocean:
    object.rpos = world_space_pos(Ocean_Frame,object.pos)

    if object.rpos.z < 0:
        object.beta = diff_angle(vector(1,0,0),vector(object.rpos.x,0,object.rpos.z))
    if object.rpos.z <= 0:
        object.beta = -diff_angle(vector(1,0,0),vector(object.rpos.x,0,object.rpos.z))

    if object.beta > 0 :
        High_Tide_Grav = High_Tide_Grav + GMm/mag2(Moon.rposobject.rpos)
    if object.beta < 0 and object.beta > -pi/2:
        ow_Tide_Grav = Low_Tide_Grav + GMm/mag2(Moon.rpos-object.rpos)

Torque_Ratio = High_Tide_Grav/Low_Tide_Grav
    #---------------------------------------------------------------------

#Computes the ratio of molecules in high tide and low tide
High_Tide = 0
Low_Tide = 0
for object in Ocean:
    if object.rpos.z < 0 and object.rpos.x > 0:
        High_Tide = High_Tide + 1
    if object.rpos.z > 0 and object.rpos.x > 0:
        Low_Tide = Low_Tide + 1
Tide_Ratio = High_Tide/Low_Tide
#---------------------------------------------------------------------
Ocean_Water = 6e19 #kg of water in the ocean
Tide_Change = 1.3 #average tide change in meters
Ocean_Depth = 1000 #average depth of ocean
Tide_Water = Ocean_Water*Tide_Change/2/Ocean_Depth
Quad_1_Tide_Mass = (Tide_Water/2)*(High_Tide/(High_Tide+Low_Tide))
Quad_4_Tide_Mass = (Tide_Water/2)*(Low_Tide/(High_Tide+Low_Tide))
#---------------------------------------------------------------------
for i in arange(0,10000,10):
    High_Tide_Point = vector((Earth_Radius-Ocean_Depth/2)*cos(pi/
    (50+i)),0,-(Earth_Radius-Ocean_Depth/2)*sin(pi/(50+i)))
    Low_Tide_Point = vector((Earth_Radius-Ocean_Depth/2)*cos(pi/
    (50+i)),0,(Earth_Radius-Ocean_Depth/2)*sin(pi/(50+i)))
    Quad_1_Pull = G*Quad_1_Tide_Mass*Moon_Mass/mag2(Moon.posHigh_Tide_Point)*norm(vector(384400,0,0)-High_Tide_Point)
    Quad_4_Pull = G*Quad_4_Tide_Mass*Moon_Mass/mag2(Moon.posLow_Tide_Point)*norm(vector(384400,0,0)-Low_Tide_Point)
    Tidal_Force = Quad_1_Pull.z + Quad_4_Pull.z
    Tidal_Torque = Tidal_Force*(Earth_Radius-Ocean_Depth/2)
    cycles = Earth.p/Tidal_Torque 
    seconds = cycles*dt
    years = seconds/(60*60*24*365)
    K.plot(pos=(pi/(50+i),years))
#---------------------------------------------------------------------
