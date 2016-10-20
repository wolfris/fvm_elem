#!/usr/bin/python -u


import numpy as np
import scipy as sp
import math



class Nodes:
    def __init__(self, nodeId, x_coord, y_coord):
        self.x = x_coord
        self.y = y_coord
        self.mid_x = self.x/2
        self.mid_y = self.y/2
        self.n_id = nodeId
        self.norm = math.sqrt(self.x**2+self.y**2)

    def diff_vec(self,r):
        vec = Nodes(0,0,0)
        vec.x = self.x - r.x
        vec.y = self.y - r.y
        vec.mid_x = r.x+0.5*(self.x - r.x)
        vec.mid_y = r.y+0.5*(self.y - r.y)
        return(vec)
        
    def normal_vec(self):
        vec = Nodes(0,0,0)
        vec.x = -self.y
        vec.y = self.x
        return(vec)

    def calc_norm(self):
        self.norm = math.sqrt(self.x**2+self.y**2)

    def dot_prod(self,r):
        a = self.x*r.x + self.y*r.y
        return(a)

    def invert_vec(self):
        self.x = -self.x
        self.y = -self.y

    def normalise(self):
        self.calc_norm()
        self.x = self.x/self.norm
        self.y = self.y/self.norm

    def is_equal(self,r,eps):
        if abs(self.x-r.x) < eps and abs(self.y-r.y) < eps:
            return 1
        else:
            return 0




class Elements:
    def __init__(self, elementId, eletype, flag, bdid, n_nodes, nodes):
        self.e_id = elementId
        self.elementType = eletype
        self.neighbours = []
        self.dx = []
        self.dy = []
        self.ds = []
        self.grad_dir = []
        self.lam = []
        self.edges = []
        self.norms = []
        self.edge_nodes = []
        self.MidpointToCentroid = []
        self.grad_dir = []
        self.boundaryId = bdid
        self.boundaryFlux = 0.0
        self.n_nodes = n_nodes
        self.nodes = nodes
        self.area = 0.0

        c_x = 0
        c_y = 0
        for i in range(self.n_nodes):
            c_x += self.nodes[i].x
            c_y += self.nodes[i].y
        c_x = c_x/self.n_nodes
        c_y = c_y/self.n_nodes
        self.centroid = Nodes(1,c_x,c_y)
        
    def find_neighbours(self,elem):
        counter = 0
        nds=[]
        nds2=[]
        for i in range(self.n_nodes):
            for j in range(elem.n_nodes):
                if (self.nodes[i].n_id == elem.nodes[j].n_id):
                    counter += 1
                    nds.append(i)
                    nds2.append(j)
        if counter == 2:
            if elem.elementType == 2 or elem.elementType == 3:
                self.neighbours.append(elem.e_id)
                elem.neighbours.append(self.e_id)
            self.edges.append(self.nodes[nds[0]].diff_vec(self.nodes[nds[1]]))
            self.edge_nodes.append(nds)
            elem.edges.append(self.nodes[nds[0]].diff_vec(self.nodes[nds[1]]))
            elem.edge_nodes.append(nds2)


            
    def calc_dx(self,elem):
        dx = self.centroid.x - elem.centroid.x
        dy = self.centroid.y - elem.centroid.y
        self.ds.append(math.hypot(dx,dy))
        self.dx.append(dx)
        self.dy.append(dy)

    def calc_area(self):
        if self.elementType == 2:
            self.area = abs(0.5*(self.edges[0].x*(-self.edges[2].y) - self.edges[0].y*(-self.edges[2].x)))
        elif self.elementType == 3:
            aux_points = []
            aux_points.append(Nodes(1,self.nodes[self.edge_nodes[0][0]].x,self.nodes[self.edge_nodes[0][0]].y))
            aux_points.append(Nodes(2,self.nodes[self.edge_nodes[0][1]].x,self.nodes[self.edge_nodes[0][1]].y))
            aux_points.append(Nodes(3,self.nodes[self.edge_nodes[1][0]].x,self.nodes[self.edge_nodes[1][0]].y))
            aux_points.append(Nodes(4,self.nodes[self.edge_nodes[1][1]].x,self.nodes[self.edge_nodes[1][1]].y))
            sum_aux = (aux_points[0].is_equal(aux_points[2],0.00001)+aux_points[0].is_equal(aux_points[3],0.00001)+
                       aux_points[1].is_equal(aux_points[2],0.00001)+aux_points[1].is_equal(aux_points[3],0.00001))
            if sum_aux == 1: 
                self.area = (abs(0.5*(self.edges[0].x*(-self.edges[1].y) - self.edges[0].y*(-self.edges[1].x))) +
                             abs(0.5*(self.edges[2].x*(-self.edges[3].y) - self.edges[2].y*(-self.edges[3].x))))
            elif sum_aux == 0:
                self.area = (abs(0.5*(self.edges[0].x*(-self.edges[2].y) - self.edges[0].y*(-self.edges[2].x))) +
                             abs(0.5*(self.edges[1].x*(-self.edges[3].y) - self.edges[1].y*(-self.edges[3].x))))
        
            
    def calc_normals(self):
        if self.elementType == 2:
            aux = 3
        elif self.elementType == 3:
            aux = 4
        elif self.elementType == 1:
            aux = 1
        for i in range(aux):
            vec_aux1 = self.edges[i].normal_vec()
            vec_aux2 = Nodes(0,vec_aux1.x,vec_aux1.y)
            vec = Nodes(0,vec_aux2.x/vec_aux2.norm,vec_aux2.y/vec_aux2.norm)

            self.norms.append(vec)
            self.edges[i].calc_norm()
            


    def calc_edgeMidpoint_centroid(self):
        if self.elementType == 2 or self.elementType == 3:
            for i in range(len(self.edges)):
                vec_aux = Nodes(0,self.edges[i].mid_x,self.edges[i].mid_y)
                vec = self.centroid.diff_vec(vec_aux)
                self.MidpointToCentroid.append(vec)
                self.MidpointToCentroid[i].calc_norm()

    def calc_lambda(self,elem,ind):
        aux_node1 = self.nodes[self.edge_nodes[ind][0]]
        aux_node2 = self.nodes[self.edge_nodes[ind][1]]
        cent_dist1 = aux_node1.diff_vec(self.centroid)
        cent_dist2 = aux_node2.diff_vec(self.centroid)
        cent_dist1.calc_norm()
        cent_dist2.calc_norm()
        semiper = (cent_dist1.norm + cent_dist2.norm + self.edges[ind].norm)/2
        x_e  = 2*math.sqrt(abs(semiper*(semiper-cent_dist1.norm)*(semiper-cent_dist2.norm)*(semiper-self.edges[ind].norm)))/self.edges[ind].norm
        x_Ex_P = self.centroid.diff_vec(elem.centroid)
        x_Ex_P.calc_norm()
        self.lam.append(x_e/x_Ex_P.norm)
