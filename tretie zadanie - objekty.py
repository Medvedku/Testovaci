import numpy as np

np.set_printoptions(precision=4)

spiel = True

class Node:
    num_of_nodes = 0
    def __init__(self, xz: list, alfa_ = 0):
        self.nd_id    = Node.num_of_nodes
        self.u        = Node.num_of_nodes*2 + 0
        self.v        = Node.num_of_nodes*2 + 1
        self.alfa     = alfa_
        self.coor     = xz

        Node.num_of_nodes += 1

class Element:
    num_of_elements = 0
    def __init__(self, start_, end_, A, E = 210*1e6):
        self.el_id       = Element.num_of_elements
        self.start       = start_.nd_id
        self.end         = end_.nd_id
        self.start_coor  = start_.coor
        self.end_coor    = end_.coor
        self.E           = E
        self.A           = A

        self.c_n = [start_.u, start_.v, end_.u, end_.v]

        self.Lx = self.end_coor[0] - self.start_coor[0]
        self.Ly = self.end_coor[1] - self.start_coor[1]
        self.L = ( (self.Lx**2) + (self.Ly**2) )**(1/2)
        assert (self.L!=0), "Element s nulovou dlzkou"

        if  self.Lx > 0  and self.Ly > 0:
            self.alfa =  np.arctan( (self.Ly/self.Lx))
        elif self.Lx < 0 and self.Ly > 0:
            self.alfa =  np.pi - np.arctan( (abs(self.Ly)/abs(self.Lx)) )
        elif self.Lx < 0 and self.Ly < 0:
            self.alfa = np.pi + np.arctan( (abs(self.Ly)/abs(self.Lx)) )
        elif self.Lx > 0 and self.Ly < 0:
            self.alfa = 2*np.pi - np.arctan(abs(self.Ly)/abs(self.Lx))
        elif self.Lx > 0 and self.Ly == 0:
            self.alfa = 0*np.pi
        elif self.Lx < 0 and self.Ly == 0:
            self.alfa = 1*np.pi
        elif self.Lx == 0 and self.Ly > 0:
            self.alfa = 0.5*np.pi
        elif self.Lx == 0 and self.Ly < 0:
            self.alfa = 1.5*np.pi
        else:
            print("Daco je velmi zle")
            quit()

        self.alfa_i = self.alfa - start_.alfa
        self.alfa_j = self.alfa -   end_.alfa

        Element.num_of_elements += 1

        self.k_elem = ((self.E * self.A)/self.L) * np.array(
        [   [ + (np.cos(self.alfa_i))*(np.cos(self.alfa_i)) , + (np.cos(self.alfa_i))*(np.sin(self.alfa_i)) , - (np.cos(self.alfa_i))*(np.cos(self.alfa_j)) , - (np.sin(self.alfa_j))*(np.cos(self.alfa_i)) ],
            [ + (np.cos(self.alfa_i))*(np.sin(self.alfa_i)) , + (np.sin(self.alfa_i))*(np.sin(self.alfa_i)) , - (np.sin(self.alfa_i))*(np.cos(self.alfa_j)) , - (np.sin(self.alfa_i))*(np.sin(self.alfa_j)) ],
            [ - (np.cos(self.alfa_j))*(np.cos(self.alfa_i)) , -((np.cos(self.alfa_j))*(np.sin(self.alfa_i))), + (np.cos(self.alfa_j))*(np.cos(self.alfa_j)) , + (np.cos(self.alfa_j))*(np.sin(self.alfa_j)) ],
            [ -((np.sin(self.alfa_j))*(np.cos(self.alfa_i))), - (np.sin(self.alfa_j))*(np.sin(self.alfa_i)) , + (np.cos(self.alfa_j))*(np.sin(self.alfa_j)) , + (np.sin(self.alfa_j))*(np.sin(self.alfa_j)) ] ]
            )

        self.T_elem = np.array(
        [   [ (np.cos(self.alfa_i)), (np.sin(self.alfa_i)),                     0,                     0 ],
            [                     0,                     0, (np.cos(self.alfa_j)), (np.sin(self.alfa_j)) ] ] )


        self.K_el = ((self.E * self.A)/self.L) * np.array(
                    [   [ +1, -1 ],
                        [ -1, +1 ] ] )

    def app_deformation(self, defs_: np.array, multiplier = 10):
        self.def_coor     = [self.start_coor[0], self.start_coor[1],self.end_coor[0], self.end_coor[1] ]
        for i in range(len(self.def_coor)):
            self.def_coor[i] += defs_[i] * multiplier

# >>>>>>>>>> Vstupy <<<<<<<<<<<<

Nodes = []
Nodes.append(Node([ 0.0, 0.0 ]))                         # 0
Nodes.append(Node([ 2.5, 0.0 ]))                         # 1
Nodes.append(Node([ 5.0, 0.0 ]))                         # 2
Nodes.append(Node([ 0.0, 1.0 ]))                         # 3
Nodes.append(Node([ 2.5, 1.5 ], np.deg2rad(230)))        # 4
Nodes.append(Node([ 5.0, 1.5 ]))                         # 5

Elements = []
Elements.append(Element( Nodes[0], Nodes[1], 3600/1e6))
Elements.append(Element( Nodes[1], Nodes[2], 3600/1e6))
Elements.append(Element( Nodes[0], Nodes[3], 2500/1e6))
Elements.append(Element( Nodes[1], Nodes[3], 2500/1e6))
Elements.append(Element( Nodes[1], Nodes[4], 2500/1e6))
Elements.append(Element( Nodes[5], Nodes[1], 2500/1e6))
Elements.append(Element( Nodes[4], Nodes[2], 2500/1e6))
Elements.append(Element( Nodes[5], Nodes[2], 2500/1e6))
Elements.append(Element( Nodes[3], Nodes[4], 3600/1e6))
Elements.append(Element( Nodes[4], Nodes[5], 3600/1e6))


# >>>>>>>>>> Začiatok Výpočtu <<<<<<<<<<<<

n = Elements[-1].el_id+1

K = np.zeros( (Nodes[-1].v+1, Nodes[-1].v+1) )
for i in Elements:
    K_temp = np.zeros( (Nodes[-1].v+1, Nodes[-1].v+1) )
    for j in range(4):
        for k in range(4):
            K_temp[i.c_n[j],i.c_n[k]] = i.k_elem[j,k]
    K += K_temp

load    = [ 0 for i in range(Nodes[-1].v+1)]

load[3] = -10
load[6] = -22 * np.cos(np.deg2rad(60))
load[7] = -22 * np.sin(np.deg2rad(60))

indexes = list(range(Nodes[-1].v+1))

deleto = [ 0, 1, 4, 9, 10, 11]
K       = np.delete(K, deleto, axis = 0)
K       = np.delete(K, deleto, axis = 1)
load    = np.delete(load, deleto, axis = 0)
indexes = np.delete(indexes, deleto, axis = 0)


delta = np.linalg.inv(K)

r_tot = np.matmul( delta, load )

r_eg = []
for i in Elements:
    r_e = []
    for j in i.c_n:
        if j in list(indexes):
            coor_ = list(indexes).index(j)
            r_e.append( r_tot[coor_] )
        else:
            r_e.append( 0 )
    r_eg.append(r_e)

r_el = [ np.matmul(Elements[i].T_elem, r_eg[i]) for i in range(n) ]

R_el = [ np.matmul(Elements[i].K_el, r_el[i]) for i in range(n) ]

S = np.array([0, 1])

N = [ np.matmul(S, R_el[i]) for i in range(n) ]

if spiel:

    for i in range(n):
        print("Osová sila v prúte {} je {} kN.".format(format(i,".3g"), format(N[i],".5g") ))

if spiel:
    import matplotlib.pyplot as plt

    for i in range(n):
        Elements[i].app_deformation(r_eg[i], 2500)

    for i in Elements:
        xs_ = np.array([i.start_coor[0], i.end_coor[0]])
        zs_ = np.array([i.start_coor[1], i.end_coor[1]])
        plt.plot(xs_, zs_, c="b", marker = "o")

    for i in Elements:
        xs = np.array([i.def_coor[0], i.def_coor[2]])
        zs = np.array([i.def_coor[1], i.def_coor[3]])
        plt.plot(xs, zs, c="r", marker = "o")

    plt.show()



# if spiel:
#     for i in Elements:
#         print("El", i.el_id,
#         " Zac. a kon. uzol",i.start, i.end,
#         " Sur", i.start_coor, i.end_coor,
#         " Dlzka", format(i.L, ".4g"),
#         " Natoc", format(np.rad2deg(i.alfa), ".6g"),
#         " alfa_i", format(np.rad2deg(i.alfa_i), ".6g"),
#         " alfa_j", format(np.rad2deg(i.alfa_j), ".6g")
#         )
#         print("kodové čísla:", i.c_n)
        # print(i.k_elem)
