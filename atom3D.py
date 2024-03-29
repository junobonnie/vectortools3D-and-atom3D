# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:59:17 2024

@author: replica
"""

import math as m
from vectortools3D import *
import pygame as pg
import sys
import itertools as it
import h5py

edges = []
for i in it.combinations(it.product([0,1],repeat=3),2):
    if sum([abs(i[0][j]-i[1][j]) for j in range(3)])==1:
        edges.append(i)
             
class Element:
    def __init__(self, name, mass, radius, color):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.color = color

    def __str__(self):
        return ('Element(name = ' + self.name + ', mass = ' + str(self.mass) +
                ', radius = ' + str(self.radius) +
                ', color = ' + str(self.color) + ')')

class Atom:
    def __init__(self, element, pos, vel = Vector(0, 0, 0)):
        self.element = element
        self.pos = pos
        self.vel = vel

    def __str__(self):
        return 'Atom(element = ' + self.element.name + ', pos(' + str(self.pos.x) + ', ' + str(self.pos.y) + ', ' + str(self.pos.z) +'), vel(' + str(self.vel.x) + ', ' + str(self.vel.y) + ', ' + str(self.vel.z) + '))'

    def is_collision(self, other):
        if isinstance(other, Atom):
            if not self == other:
                d = self.pos - other.pos
                v1 = self.vel
                v2 = other.vel

                v1_ = v1.dot(d)*d/(d.dot(d))
                v2_ = v2.dot(d)*d/(d.dot(d))

                return (d.dot(d) < (self.element.radius + other.element.radius)**2) and (d.dot(v1_-v2_) < 0)
            else:
                return False

    def collision(self, other, dt):
        if isinstance(other, Atom):
            if not self == other:
                d = self.pos - other.pos
                m1 = self.element.mass
                v1 = self.vel
                m2 = other.element.mass
                v2 = other.vel

                v1_ = v1.dot(d)*d/(d.dot(d))
                v2_ = v2.dot(d)*d/(d.dot(d))

                if (d.dot(d) < (self.element.radius + other.element.radius)**2) and (d.dot(v1_-v2_) < 0):
                    self.pos -= self.vel*dt
                    other.pos -= other.vel*dt
                    v1__ = (m1-m2)/(m1+m2)*v1_ + 2*m2/(m1+m2)*v2_
                    v2__ = 2*m1/(m1+m2)*v1_ + (m2-m1)/(m1+m2)*v2_
                    self.vel = v1 - v1_ + v1__
                    other.vel = v2 - v2_ + v2__

    def fusion(self, other_Atom):
        new_Atom = None
        if not self == other_Atom:
            d = self.pos - other_Atom.pos
            if (d.dot(d) < (self.element.radius + other_Atom.element.radius)**2):
                new_element = Element(name = 'New atom', mass = self.element.mass + other_Atom.element.mass, 
                                      radius = m.sqrt(self.element.radius**2 + other_Atom.element.radius**2),
                                      color = self.element.color + other_Atom.element.color)
                new_Atom = Atom(element = new_element, 
                                pos = (self.element.mass*self.pos + other_Atom.element.mass*other_Atom.pos)/(self.element.mass + other_Atom.element.mass),
                                vel = (self.element.mass*self.vel + other_Atom.element.mass*other_Atom.vel)/(self.element.mass + other_Atom.element.mass))
        return new_Atom
    
class World:
    def __init__(self, t, atoms, gravity):
        self.t = t
        self.atoms = atoms
        self.gravity = gravity

    def __str__(self):
        return ('World(t = ' + str(self.t) + ', atoms = ' + str(self.atoms) 
                + ', gravity = ' + str(self.gravity) + ')')

class Render:
    def __init__(self, screen, width, height, depth, angle=(0,0,0), omega = (0,0,0), focus_factor = 1):
        if not screen == None:
            pg.init()
        self.screen = screen
        self.width = width
        self.height = height
        self.depth = depth
        self.focus_factor = focus_factor
        self.render_vector = Vector(0, height, 0)
        self.render_metric = Tensor(1, 0, 0, 
                                    0, -1, 0,
                                    0, 0, 1)
        self.origin_vector = Vector(width/2, height/2, depth/2)
        self.angle = angle
        self.omega = omega
        self.t = 0

    def update_time(self, t):
        self.t = t

    def render_SO3(self):
        if self.omega == (0, 0, 0):
            omega_rotaion = Tensor(1, 0, 0,
                                    0, 1, 0,
                                    0, 0, 1)
        else:
            omega_rotaion = SO3_x(self.omega[0]*self.t).dot(
                            SO3_y(self.omega[1]*self.t).dot(
                            SO3_z(self.omega[2]*self.t)))
        if self.angle == (0, 0, 0):
            angle_rotaion = Tensor(1, 0, 0,
                                    0, 1, 0,
                                    0, 0, 1)
        else:
            angle_rotaion = SO3_x(self.angle[0]).dot(
                            SO3_y(self.angle[1]).dot(
                            SO3_z(self.angle[2])))
        return omega_rotaion.dot(angle_rotaion)

    def rendering_vector(self, vector):
        return self.render_vector + self.render_metric.dot(vector + self.origin_vector)

    def text(self, text, font, size, pos, color):
        pos = self.render_SO3().dot(pos)
        value = self.focus_factor*self.depth/(pos.z + self.depth)
        font_ = pg.font.SysFont(font, int(value*size))
        text_ = font_.render(text, True, color)
        render_pos = self.rendering_vector(value*pos)
        self.screen.blit(text_, (render_pos.x, render_pos.y))

    def polygon(self, positions, color):
        P_list = []
        for pos in positions:
            pos = self.render_SO3().dot(pos)
            value = self.focus_factor*self.depth/(pos.z + self.depth)
            P = self.rendering_vector(value*pos)
            P_list.append([P.x, P.y])
        pg.draw.aalines(self.screen, color, True, P_list, True)
        
    def cube(self, pos, width, height, depth, color):
        for ee in edges:
            edge = []
            for e in ee:
                edge.append(pos + Vector(e[0]*width, e[1]*height, e[2]*depth))
            self.polygon(edge, color)
            
    def circle(self, pos, radius, color):
        pos = self.render_SO3().dot(pos)
        value = self.focus_factor*self.depth/(pos.z + self.depth)
        render_pos =self.rendering_vector(value*pos)
        pg.draw.circle(self.screen, color, (render_pos.x, render_pos.y), value*radius)

    def atom(self, atom):
        self.circle(atom.pos, atom.element.radius, atom.element.color)

class Simulator:
    def __init__(self, dt, world, render, grid_size = 100):
        self.dt = dt
        self.world = world
        self.render = render
        self.count_screen = 0
        self.count_snapshot = 0
        self.region_size = 0
        self.grid_size = grid_size
        self.grid = None

    def clock(self):
        self.world.t = self.world.t + self.dt
        return self.world.t

    def draw_background(self, color):
        self.render.screen.fill(color)

    def draw_grid(self, unit_size):
        grey = (200, 200, 200)
        for x in range(0, int(self.render.width/2), unit_size):
            for z in range(0, int(self.render.depth/2), unit_size):
                for sign in ((-1,-1),(1,-1),(-1,1),(1,1)):
                    self.render.polygon([Vector(sign[0]*x, -self.render.height/2, sign[1]*z), Vector(sign[0]*x, self.render.height/2, sign[1]*z)], grey)
        for x in range(0, int(self.render.width/2), unit_size):
              for y in range(0, int(self.render.height/2), unit_size):
                  for sign in ((-1,-1),(1,-1),(-1,1),(1,1)):
                      self.render.polygon([Vector(sign[0]*x, sign[1]*y, -self.render.depth/2), Vector(sign[0]*x, sign[1]*y, self.render.depth/2)], grey)
        for y in range(0, int(self.render.height/2), unit_size):
              for z in range(0, int(self.render.depth/2), unit_size):
                  for sign in ((-1,-1),(1,-1),(-1,1),(1,1)):
                      self.render.polygon([Vector(-self.render.width/2, sign[0]*y, sign[1]*z), Vector(self.render.width/2, sign[0]*y, sign[1]*z)], grey)
                
    def draw_atom(self):
        for atom in self.world.atoms:
            if ((-self.render.width/2 < atom.pos.x < self.render.width/2) and 
                (-self.render.height/2 < atom.pos.y < self.render.height/2) and
                (-self.render.depth/2 < atom.pos.z < self.render.depth/2)):
                self.render.atom(atom)

    def update_region_size(self):
        self.region_size = 0
        for atom in self.world.atoms:
            for p in atom.pos.list():
                self.region_size = max(self.region_size, 2*abs(p))
        return self.region_size

    def make_grid(self):
        n = int(self.region_size//self.grid_size+1)
        grid = [[] for i in range(n**3)]
        for atom in self.world.atoms:
            i = int((self.region_size/2 + atom.pos.x)//self.grid_size)
            j = int((self.region_size/2 + atom.pos.y)//self.grid_size)
            k = int((self.region_size/2 + atom.pos.z)//self.grid_size)
            if (0 <= i < n) and (0 <= j < n) and (0 <= k < n):
                grid[i+n*j+n*n*k].append(atom)
        self.grid = grid
    
    def get_near_atoms(self, atom):
        n = int(self.region_size//self.grid_size+1)
        i = int((self.region_size/2 + atom.pos.x)//self.grid_size)
        j = int((self.region_size/2 + atom.pos.y)//self.grid_size)
        k = int((self.region_size/2 + atom.pos.z)//self.grid_size)
        atoms = []
        for i_ in (i-1, i, i+1):
            for j_ in (j-1, j, j+1):
                for k_ in (k-1, k, k+1):
                    if (0 <= i_ < n) and (0 <= j_ < n) and (0 <= k_ < n):
                        atoms += self.grid[i_+n*j_+n*n*k_]
        return atoms

    def atom_atom_collision(self):
        self.update_region_size()
        self.make_grid()
        for atom in self.world.atoms:
            atoms = self.get_near_atoms(atom)
            for other_atom in atoms:
                #self.render.polygon([atom.pos, other_atom.pos], red)
                atom.collision(other_atom, self.dt)
    
    def atom_atom_fusion(self):
        while True:
            for atom in self.world.atoms[:]:
                for other_atom in self.world.atoms[:]:
                    new_atom = atom.fusion(other_atom)
                    if not new_atom == None:
                        self.world.atoms.remove(atom)
                        self.world.atoms.remove(other_atom)
                        self.world.atoms.append(new_atom)
                        break
                if not new_atom == None:
                    break
            if new_atom == None:
                break
                
    def main(self):
        x_ = []
        v_ = []
        for atom in self.world.atoms:
            new_v = atom.vel + self.world.gravity*self.dt
            v_.append(new_v)
            x_.append(atom.pos + new_v*self.dt)

        count = 0
        for atom in self.world.atoms:
            atom.pos = x_[count]
            atom.vel = v_[count]
            count = count + 1

    def save_screen(self, directory, skip_number = 0):
        if self.count_screen%(skip_number+1) == 0:
            img = directory + '/%08d.png' % (self.count_screen)
            pg.image.save(self.render.screen, img)
        self.count_screen += 1
        
    # def save_snapshot(self, directory, skip_number = 0):
    #     if self.count_snapshot%(skip_number+1) == 0:
    #         snapshot = directory + '/snapshot_%08d.txt' % (self.count_snapshot)
    #         with open(snapshot, "w") as f:
    #             atoms_info = ''
    #             count = 0
    #             for atom in self.world.atoms:
    #                 count += 1
    #                 atoms_info = atoms_info + 'atom' + str(count) + '{ element{ name:' + atom.element.name + ', mass:' + str(atom.element.mass) + ', radius:' + str(atom.element.radius) + ', color:' + str(atom.element.color) + ' }' + ', pos:' + str(atom.pos) + ', vel:' + str(atom.vel) + ' }, '
                    
    #             f.write('world{ t:' + str(self.world.t) + ', gravity:' + str(self.world.gravity) + ', atoms{ ' + atoms_info + ' }' + ' }')
    #     self.count_snapshot += 1

    def save_snapshot(self, directory, skip_number = 0):
        if self.count_snapshot%(skip_number+1) == 0:
            snapshot = directory + '/snapshot_%08d.hdf5' % (self.count_snapshot)
            f = h5py.File(snapshot, 'w')
            world = f.create_group('world')
            world.attrs['t'] = self.world.t
            world.attrs['gravity'] = self.world.gravity.list()
            atoms = world.create_group('atoms')
            count = 0
            for atom in self.world.atoms:
                count += 1
                atom_ = atoms.create_group('atom'+str(count))
                atom_.attrs['element'] = atom.element.name
                atom_.attrs['mass'] = atom.element.mass
                atom_.attrs['radius'] = atom.element.radius
                atom_.attrs['color'] = atom.element.color
                atom_.attrs['pos'] = atom.pos.list()
                atom_.attrs['vel'] = atom.vel.list()
            f.close()
        self.count_snapshot += 1
        
    # def load_snapshot(self, snapshot_file):
    #     with open(snapshot_file, "r") as f:
    #         snapshot = f.read()
    #         snapshot = snapshot.replace('world{ t:', '').replace(', gravity:', '#').replace(', atoms', '#').replace(' },  } }', '') 
    #         snapshot = snapshot.split('#')
    #         t = float(snapshot[0])
    #         gravity = eval(snapshot[1])
            
    #         atoms_raw = snapshot[2]
    #         atoms_raw = atoms_raw.replace('{ atom', '').split('atom')
    #         atoms = []
    #         for atom in atoms_raw:
    #             atom = atom.replace('element{ ', '#').replace(' }, pos:', '#').replace(', vel:', '#').replace(' }, ', '')
    #             atom = atom.split('#')
    #             try:
    #                 element_raw = atom[1]
    #                 element_raw = element_raw.replace('name:', '#').replace(', mass:', '#').replace(', radius:', '#').replace(', color:', '#')
    #                 element_raw = element_raw.split('#')
    #                 atoms.append(Atom(Element(element_raw[1], float(element_raw[2]), float(element_raw[3]), eval(element_raw[4])), eval(atom[2]), eval(atom[3])))
    #             except:
    #                 pass
    #     self.world = World(t, atoms, gravity)
        
    def list_to_vector(self, list):
        return Vector(list[0], list[1], list[2])

    def load_snapshot(self, snapshot_file):
        f = h5py.File(snapshot_file, 'r')
        world = f['world']
        t = world.attrs['t']
        gravity = self.list_to_vector(world.attrs['gravity'])
        atoms = []
        for atom_ in world['atoms']:
            atom = world['atoms'][atom_]
            element = Element(atom.attrs['element'], atom.attrs['mass'], atom.attrs['radius'], atom.attrs['color'])
            pos = self.list_to_vector(atom.attrs['pos'])
            vel = self.list_to_vector(atom.attrs['vel'])
            atoms.append(Atom(element, pos, vel))
        self.world = World(t, atoms, gravity)
        f.close()
                
if __name__ == '__main__':
    width = 1000
    height = 1000
    depth = 1000

    screen = pg.display.set_mode((width, height))
    render = Render(screen, width, height, depth, (m.pi/4, m.pi/4, 0), (0, 10, 0))
    #clock = pg.time.Clock()

    black = pg.Color('black')
    white = pg.Color('white')
    red = pg.Color('red')
    green = pg.Color('green')
    blue = pg.Color('blue')

    e1 = Element(name = 'Helium', mass = 1, radius = 10, color = red)
    atom1 = Atom(e1, Vector(-200, 0, 1), Vector(500, 0, 0))
    atom2 = Atom(e1, Vector(0, 0, 0))
    atom3 = Atom(e1, Vector(25, -10, 0))
    atom4 = Atom(e1, Vector(25, 10, 0))
    atom5 = Atom(e1, Vector(50, -20, 0))
    atom6 = Atom(e1, Vector(50, 0, 0))
    atom7 = Atom(e1, Vector(50, 20, 0))
    
    atoms = [atom1, atom2, atom3, atom4, atom5, atom6, atom7]

    gravity = Vector(0, -10, 0)*0

    world = World(0, atoms, gravity)

    simulator = Simulator(0.01, world, render)
    simulator.load_snapshot('snapshots/snapshot_00000000.hdf5')
    while True:
        t = simulator.clock()
        simulator.render.update_time(t)
        simulator.draw_background(white)
        simulator.draw_grid(100)
        simulator.atom_atom_collision()
        #simulator.atom_atom_fusion()
        simulator.main()
        simulator.draw_atom()
        
        render.cube(Vector(-100, -250 ,0), 200, 500, 100, blue)
        render.text('pos = (%.2f, %.2f, %.2f)'%(atom1.pos.x, atom1.pos.y, atom1.pos.z) , None, 30, Vector(atom1.pos.x -100, atom1.pos.y - 30, atom1.pos.z), black)
        render.text('vel = (%.2f, %.2f, %.2f)'%(atom1.vel.x, atom1.vel.y, atom1.vel.z) , None, 30, Vector(atom1.pos.x -100, atom1.pos.y - 50, atom1.pos.z), black)

        render.text('pos = (%.2f, %.2f, %.2f)'%(atom7.pos.x, atom7.pos.y, atom7.pos.z) , None, 30, Vector(atom7.pos.x -100, atom7.pos.y - 30, atom1.pos.z), blue)
        render.text('vel = (%.2f, %.2f, %.2f)'%(atom7.vel.x, atom7.vel.y, atom7.vel.z) , None, 30, Vector(atom7.pos.x -100, atom7.pos.y - 50, atom1.pos.z), blue)

        for event in pg.event.get():
            if event.type == pg.QUIT:
                sys.exit()
        #clock.tick(100)
        pg.display.update()
        
        #simulator.save_screen('images/pocket_ball_demo')
        simulator.save_snapshot('snapshots', 99)
