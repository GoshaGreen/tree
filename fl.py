from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import svg # pip install svg.py
import os


NORM = 700 / (872.3-0)
WOOD_WIDTH = 18
PLATE_SIZE = 35

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

class Point():
    def __init__(self, name, posX, posY):
        self.name = name
        self.x = posX
        self.y = posY
        self.xy = np.array([self.x, self.y])
    def xy(xy, name =''):
        return Point(name, xy[0], xy[1])
    def distance(self, other)->int:
        if not isinstance(other, Point): raise
        return np.linalg.norm(other.xy-self.xy)
    def __str__(self):
        return f'{self.name}[{self.x}:{self.y}]'
    def __repr__(self):
        return self.__str__()

# Tree Node
class TN:
    IND = ' |'
    # Connection types
    C_STOP = 'C_STOP'
    C_DIRECT = 'C_DIRECT'
    C_LICK = 'C_LICK'
    C_NORMAL = 'C_NORMAL' # should be the first (0st) child 
    C_FORCE = 'C_FORCE'
    C_STOP_LICK = 'C_STOP_LICK'
    C_FLOW = 'C_FLOW'
    C_GLUE = 'C_GLUE'

    def __init__(self, point: Point, connType = C_STOP):
        self.name = point.name
        self.parent = None
        self.point = point
        self.connType = connType
        self.childs = list()
    def addChild(self, connType, child):
        self.childs.append(child)
        self.childs[-1].setParent(self)
        self.childs[-1].connType = connType
        return self
    def setParent(self, parent):
        self.parent = parent
    def genStr(self, indient):
        stringos = ''
        clds = len(self.childs)
        stringos += f'{indient}`{self.name}    -----> [childs:{clds}]\n'
        for child in self.childs:
            clds -= 1
            indt = indient + TN.IND
            if clds <=0:
                indt = indt[:-2] + '  '
            stringos += child.genStr(indt) 
        return stringos
    def distance(self, other)->int:
        if isinstance(other, TN):
            return np.linalg.norm(other.point.xy-self.point.xy)
        elif isinstance(other, Point):
            return np.linalg.norm(other.xy-self.point.xy)
        else:
            print(isinstance(other, TN), isinstance(other, Point), type(Point), type(TN))
            raise

    def calcBranches(self, parBr = None):
        branches = list()
        for child in self.childs:
            branches.append(
                TreeBranch(f'{self.name:_>3}->{child.name:_<3}', self, child, parBr=parBr)
            )
            branches += child.calcBranches(branches[-1])
        return branches
    def __str__(self):
        return self.genStr('')
    def __repr__(self):
        return self.__str__()

def unit_vector(vector):
    return vector / np.linalg.norm(vector)
def angle_360(v_in):
    v = unit_vector(v_in)
    c = np.arccos(v[0])
    b = math.copysign(c, v[1])
    return np.rad2deg((0 if b > 0 else 2*np.pi) + b)
def angle_between(v1, v2):
    return angle_360(v2) - angle_360(v1)

class TreeBranch:
    inner_id = 0
    def __init__(self, name: str , node1: TN, node2: TN, parBr = None):
        self.name = name
        self.node1 = TN(node1.point, node2.connType)
        self.node2 = TN(node2.point, TN.C_STOP)
        self.node2.childs = [node2.childs[0]] if len(node2.childs) else []
        self.parBr = parBr
        self.childs = list()
        self.p11 = None
        self.p12 = None
        self.p21 = None
        self.p22 = None
        self.leng1 = 0
        self.leng2 = 0
        self.angle1 = 90
        self.angle2 = 90
        self.platePoints = list()
        self.inner_id = TreeBranch.inner_id
        TreeBranch.inner_id = TreeBranch.inner_id + 1
        if isinstance(self.parBr, TreeBranch) : self.parBr.childs.append(self)
        self.calcPoints()
    def angleParent(self):
        v1 = self.parBr.node1.point.xy - self.node1.point.xy
        v2 = self.node2.point.xy - self.node1.point.xy
        return angle_between(v2,v1)
    def angleChild(self, childIdx = 0):
        v1 = self.node1.point.xy-self.node2.point.xy
        v2 = self.node2.childs[childIdx].point.xy-self.node2.point.xy
        return angle_between(v1,v2)
    def calcPlates(self): # calc plates near node 2
        plates = list()
        noPlate = (len(self.childs) < 2)
        for child in self.childs:
            if child.node1.connType in [TN.C_DIRECT, TN.C_NORMAL, TN.C_FORCE, TN.C_FLOW]:break
        else:
            noPlate = True
        if not noPlate:
            V21 = unit_vector(self.node1.point.xy - self.node2.point.xy)
            Vort = unit_vector([-V21[1], V21[0]])
            plateName = f'{self.name}'
            self.platePoints.append(Point.xy(self.node2.point.xy+V21*PLATE_SIZE + Vort * (WOOD_WIDTH / 2), plateName))
            self.platePoints.append(Point.xy(self.node2.point.xy+V21*PLATE_SIZE - Vort * (WOOD_WIDTH / 2), plateName))
            # if 
            for child in self.childs:
                if child.node1.connType in [TN.C_STOP, TN.C_GLUE]: continue
                self.platePoints.append(Point.xy(child.p12.xy))
                V21 = unit_vector(child.node2.point.xy- child.node1.point.xy )
                Vort = unit_vector([-V21[1], V21[0]])
                self.platePoints.append(Point.xy(child.node1.point.xy+V21*PLATE_SIZE + Vort * (WOOD_WIDTH / 2), plateName))
                self.platePoints.append(Point.xy(child.node1.point.xy+V21*PLATE_SIZE - Vort * (WOOD_WIDTH / 2), plateName))
                self.platePoints.append(Point.xy(child.p11.xy))
            plates.append(self.platePoints)
        for child in self.childs:
            plates += child.calcPlates()
        return plates
    def calcPoints(self):
        baseLen = self.node1.distance(self.node2)
        dy1 = 0
        dx1 = 0
        dy2 = 0
        dx2 = 0
        
        # TN.C_STOP:
        dy1 = - (WOOD_WIDTH / 2) / baseLen * (self.node2.point.x - self.node1.point.x)
        dx1 =   (WOOD_WIDTH / 2) / baseLen * (self.node2.point.y - self.node1.point.y)
        dy2 = - (WOOD_WIDTH / 2) / baseLen * (self.node1.point.x - self.node2.point.x)
        dx2 =   (WOOD_WIDTH / 2) / baseLen * (self.node1.point.y - self.node2.point.y)
        childConnType = self.node2.childs[0].connType if len(self.node2.childs) else TN.C_STOP
        connType = self.node1.connType
        # direction changed calculations
        if connType == TN.C_DIRECT:
            # find closest bw parentBranch.p21 and parentBranch.p22
            t = 'r' if self.node2.distance(self.parBr.p22) < self.node2.distance(self.parBr.p21) else 'l'
            a_p = self.parBr.p22.xy if t == 'r' else self.parBr.p21.xy
            p_p = self.node2.point.xy
            Vap = p_p - a_p
            o_p = a_p + Vap/2
            a_a = np.rad2deg(np.arcsin(WOOD_WIDTH/2/np.linalg.norm(Vap)))
            pon_a = angle_360(Vap) + 2 * a_a * (1 if t =='r' else -1)
            Von = np.linalg.norm(Vap/2) * unit_vector(np.array([np.cos(pon_a * np.pi/180.), np.sin(pon_a * np.pi/180.)]))
            n_p = o_p + Von
            Vaf = self.node1.point.xy - self.parBr.node1.point.xy
            Van = n_p-a_p
            be_a = angle_between(Van, Vaf)
            BC = WOOD_WIDTH/2/np.sin(be_a * np.pi/180.)
            Vs = unit_vector(Vaf) * BC
            # set new self.node1 point
            node1 = a_p - Vs* (1 if t =='r' else -1)
            self.node1.point = Point(self.node1.point.name, node1[0], node1[1])
            Vb1, Vb2 = self.node1.point.xy - Vs, self.node1.point.xy + Vs
            self.p11 = Point(self.name, Vb1[0], Vb1[1])
            self.p12 = Point(self.name, Vb2[0], Vb2[1])
            dx1, dy1 = None, None
            # plt.plot([n_p[0],self.node1.point.x], [n_p[1],self.node1.point.y], 'ro')
        elif connType == TN.C_FORCE:
            prevBro = None
            for bro in self.parBr.childs:
                # print(bro.inner_id, self.inner_id, bro.inner_id == self.inner_id)
                if bro.inner_id == self.inner_id:
                    break
                if bro.node1.connType not in [TN.C_GLUE]:
                    prevBro = bro
            else:
                raise
            if prevBro == None: raise

            # find closest
            t = 'r' if self.node2.distance(prevBro.p12) < self.node2.distance(prevBro.p11) else 'l'
            a_p = prevBro.p12.xy if t == 'r' else prevBro.p11.xy
            p_p = self.node2.point.xy
            Vap = p_p - a_p
            o_p = a_p + Vap/2
            a_a = np.rad2deg(np.arcsin(WOOD_WIDTH/2/np.linalg.norm(Vap)))
            pon_a = angle_360(Vap) + 2 * a_a * (1 if t =='r' else -1)
            Von = np.linalg.norm(Vap/2) * unit_vector(np.array([np.cos(pon_a * np.pi/180.), np.sin(pon_a * np.pi/180.)]))
            n_p = o_p + Von
            Vaf = self.node1.point.xy - prevBro.node2.point.xy
            Van = n_p-a_p
            be_a = angle_between(Van, Vaf)
            BC = WOOD_WIDTH/2/np.sin(be_a * np.pi/180.)
            Vs = unit_vector(Vaf) * BC
            # set new self.node1 point
            node1 = a_p - Vs* (1 if t =='r' else -1)
            self.node1.point = Point(self.node1.point.name, node1[0], node1[1])
            Vb1, Vb2 = self.node1.point.xy - Vs, self.node1.point.xy + Vs
            self.p11 = Point(self.name, Vb1[0], Vb1[1])
            self.p12 = Point(self.name, Vb2[0], Vb2[1])
            dx1, dy1 = None, None
        elif connType == TN.C_GLUE:
            prevBro = None
            for bro in self.parBr.childs:
                if bro.inner_id == self.inner_id:
                    break
                prevBro = bro
            else:
                raise
            if prevBro == None: raise
            brlen = self.node1.distance(self.node2)
            Vab = prevBro.node1.point.xy + brlen * unit_vector(prevBro.node2.point.xy - prevBro.node1.point.xy)
            self.node1 = TN(prevBro.node1.point, self.node1.connType)
            self.node2 = TN(
                Point(self.node2.point.name, Vab[0], Vab[1]),
                self.node2.connType
            )
            self.p11 = prevBro.p11
            self.p12 = prevBro.p12
            dx1, dy1 = None, None

        # node2 points calc
        # LICK first
        if childConnType == TN.C_LICK:
            dy2 = 0
            dx2 = 0
        # C_STOP_LICK first
        elif childConnType == TN.C_STOP_LICK:
            ac = self.angleChild(0)
            c = (WOOD_WIDTH / 2) / np.sin((180-ac)*np.pi/180.)
            Vm = -unit_vector(self.node2.childs[0].point.xy-self.node2.point.xy) * c
            dx2, dy2 = Vm[0], Vm[1]
        # C_STOP_LICK second
        elif connType == TN.C_STOP_LICK:
            dy2 = 0
            dx2 = 0
        elif childConnType == TN.C_STOP:
            baseLen = self.node1.distance(self.node2)
            dy2 = - (WOOD_WIDTH / 2) / baseLen * (self.node1.point.x - self.node2.point.x)
            dx2 =   (WOOD_WIDTH / 2) / baseLen * (self.node1.point.y - self.node2.point.y)
        elif childConnType == TN.C_NORMAL:
            ac = self.angleChild(0)
            # check if child has TN.C_NORMAL 
            v1 = self.node1.point.xy - self.node2.point.xy
            bb = angle_360(v1) 
            bt = abs((WOOD_WIDTH / 2) / np.sin(ac / 2 * np.pi /180.))
            aa = bb + ac / 2
            dy2 = bt * np.sin((aa)*np.pi/180.)
            dx2 = bt * np.cos((aa)*np.pi/180.)
            # print('second', self.name, [f'{x:7.2f}' for x in [bb, ac, ac/2, bt, aa, dx2, dy2]])
        elif 0:
            pass

        # node1 points calc
        # LICK first
        if childConnType == TN.C_LICK:
            dy1 = 0
            dx1 = 0
        # LICK second
        elif connType == TN.C_LICK:
            ap = self.angleParent()
            blp = self.node1.distance(self.parBr.node1)
            ol = (WOOD_WIDTH / 2) / np.sin(abs(ap) * np.pi/180.)
            parDy = self.node1.point.y - self.parBr.node1.point.y 
            parDx = self.node1.point.x - self.parBr.node1.point.x
            dy1 = ol * parDy / blp
            dx1 = ol * parDx / blp
        # C_STOP_LICK second
        elif connType== TN.C_STOP_LICK:
            dy1 = 0
            dx1 = 0
        elif connType == TN.C_NORMAL:
            ap = self.angleParent()
            v1 = self.node2.point.xy - self.node1.point.xy
            bb = angle_360(v1)
            bt = abs((WOOD_WIDTH/2) / np.sin(ap/2*np.pi/180.)) 
            aa = ap/2+bb
            dy1 = -bt * np.sin((aa)*np.pi/180.)
            dx1 = -bt * np.cos((aa)*np.pi/180.)
            # print('first ',self.name, ['{:7.2f}'.format(x) for x in [bb, ap, ap/2, bt, aa, dx1, dy1]])
        elif connType == TN.C_FLOW:
            self.node1.point = Point(self.node1.point.name, self.parBr.node2.point.x, self.parBr.node2.point.y)
            Vad = self.node2.point.xy - self.node1.point.xy
            bb = self.angleParent() - 90
            bc = (WOOD_WIDTH/2) / np.cos( bb * np.pi/180.)
            Vm = self.node1.point.xy + unit_vector(Vad) * abs(bc)
            Voa = bc * unit_vector(self.node1.point.xy-self.parBr.node1.point.xy)
            Vb1, Vb2 = Vm + Voa, Vm - Voa
            self.p11 = Point(self.name, Vb1[0], Vb1[1])
            self.p12 = Point(self.name, Vb2[0], Vb2[1])
            dx1, dy1 = None, None
        elif 0:
            pass

        if dx1 != None and dy1 != None:
            self.p11 = Point(self.name, self.node1.point.x + dx1, self.node1.point.y + dy1)
            self.p12 = Point(self.name, self.node1.point.x - dx1, self.node1.point.y - dy1)
        if dx2 != None and dy2 != None:
            self.p21 = Point(self.name, self.node2.point.x + dx2, self.node2.point.y + dy2)
            self.p22 = Point(self.name, self.node2.point.x - dx2, self.node2.point.y - dy2)
        self.calcParams()

    def calcParams(self):
        self.leng1 = self.p12.distance(self.p22)
        self.leng2 = self.p11.distance(self.p21)

        self.angle11 = abs(angle_between(self.p12.xy-self.p11.xy, self.p22.xy-self.p11.xy))
        self.angle12 = abs(angle_between(self.p11.xy-self.p12.xy, self.p21.xy-self.p12.xy))
        self.angle22 = abs(angle_between(self.p21.xy-self.p22.xy, self.p11.xy-self.p22.xy))
        self.angle21 = abs(angle_between(self.p22.xy-self.p21.xy, self.p12.xy-self.p21.xy))

        self.angle11 = self.angle11 if self.angle11 < 180 else 360 -self.angle11 
        self.angle12 = self.angle12 if self.angle12 < 180 else 360 -self.angle12 
        self.angle21 = self.angle21 if self.angle21 < 180 else 360 -self.angle21 
        self.angle22 = self.angle22 if self.angle22 < 180 else 360 -self.angle22 

    def getPoints(self):
        return [self.p11, self.p12, self.p21, self.p22]
    def __str__(self):
        return f'[{self.inner_id}]{self.name:12} [{self.node1.point}({self.node1.connType}) - {self.node2.point}]'
    def __repr__(self):
        return self.__str__()

def mirroredPoint(name, posX, posY):
    return Point(name, posX * NORM, (2509-posY) * NORM)

A = mirroredPoint("A",669,2509)
B = mirroredPoint("B",746,2163)
C = mirroredPoint("C",702,1796)
D = mirroredPoint("D",406,1499)
E = mirroredPoint("E",409,1665)
F = mirroredPoint("F",104,1467)
G = mirroredPoint("G",0,1522)
H = mirroredPoint("H",0,1304)
I = mirroredPoint("I",472,1300)
J = mirroredPoint("J",557,1650)
K = mirroredPoint("K",757,1581)
L = mirroredPoint("L",667,1461)
M = mirroredPoint("M",812,1369)
N = mirroredPoint("N",872,1353)
O = mirroredPoint("O",738,1174)
P = mirroredPoint("P",872,1058)
Q = mirroredPoint("Q",660,966)
R = mirroredPoint("R",458,934)
S = mirroredPoint("S",344,966)
T = mirroredPoint("T",207,1002)
U = mirroredPoint("U",257,945)
V = mirroredPoint("V",305,783)
W = mirroredPoint("W",94,691)
X = mirroredPoint("X",330,657)
Y = mirroredPoint("Y",751,742)
Z = mirroredPoint("Z",617,576)
A1 = mirroredPoint("A1",873,685)
A2 = mirroredPoint("A2",0,2509)
A3 = mirroredPoint("A3",872.3,2509)


ROOT = TN(A3).addChild(
    TN.C_STOP, TN(A).addChild(
        TN.C_LICK, TN(B).addChild(
            TN.C_NORMAL, TN(C).addChild(
                TN.C_NORMAL, TN(D).addChild(
                    TN.C_NORMAL, TN(F).addChild(
                        TN.C_NORMAL, TN(H).addChild(TN.C_STOP_LICK, TN(A2))
                    ).addChild(
                        TN.C_FORCE, TN(G).addChild(TN.C_STOP_LICK, TN(A2))
                    )
                ).addChild(
                    TN.C_DIRECT, TN(I)
                )
            ).addChild(
                TN.C_GLUE, TN(J).addChild(TN.C_FLOW, TN(E))
            ).addChild(
                TN.C_FORCE, TN(M).addChild(
                    TN.C_NORMAL, TN(Q).addChild(
                        TN.C_NORMAL, TN(R).addChild(
                            TN.C_NORMAL, TN(V).addChild(
                                TN.C_NORMAL, TN(W)
                            ).addChild(
                                TN.C_FORCE, TN(X)
                            )
                        ).addChild(
                            TN.C_FORCE, TN(T)
                        ).addChild(
                            TN.C_GLUE, TN(S).addChild(TN.C_FLOW, TN(U))
                        )
                    ).addChild(
                        TN.C_FORCE, TN(Y).addChild(
                            TN.C_NORMAL, TN(A1).addChild(TN.C_STOP_LICK, TN(A3))
                        ).addChild(
                            TN.C_DIRECT, TN(Z)
                        )
                    )
                ).addChild(
                    TN.C_GLUE, TN(O).addChild(
                        TN.C_FLOW, TN(P).addChild(TN.C_STOP_LICK, TN(A3))
                    )
                ).addChild(
                    TN.C_FORCE, TN(N).addChild(TN.C_STOP_LICK, TN(A3))
                )
            ).addChild(
                TN.C_GLUE, TN(K).addChild(TN.C_FLOW, TN(L))
            )
        )
    )
)

# VV = [1,0]
# NN = 12
# for i in range(-NN,NN):
#     an1 = 2*np.pi/NN*i
#     print('{:5.2f} -> {:5.2f}'.format(np.rad2deg(an1), angle_between([1,0], [np.cos(an1),np.sin(an1)])))

print(ROOT)

branches = ROOT.calcBranches()
[print(br) for br in branches]

plates = branches[0].calcPlates()

# print()
# print(f'Branch name |  length1 |  length2 | angle11| angle12| angle21| angle22 ')
# print(f'            .          .          .        .        .        .        .         ')
# summleng = 0
# numBranchs = 0
# for br in branches:
#     if math.isnan(br.angle11) or math.isnan(br.angle12) or math.isnan(br.angle21) or math.isnan(br.angle22):
#         continue
#         pass
#     print(f'{br.name:12}.{max(br.leng1, br.leng2):8.2f}mm.{min(br.leng1, br.leng2):8.2f}mm.'+''.join([f'{a:7.1f}{'!' if not 45<a<135 else ' '}.' for a in [br.angle11,br.angle12,br.angle21,br.angle22, br.angle11+br.angle12,br.angle21+br.angle22]]))
#     summleng += max(br.leng1, br.leng2)
#     numBranchs += 1
# print()
# print(f'total length: {summleng:8.2f}mm, {numBranchs} branches,  with 10mm gap total: {summleng + 10 * numBranchs:8.2f}mm')
# print()

# print('='*80)
dname = os.path.join(SCRIPT_DIR, 'svgs') 
if not os.path.exists(dname):
    os.makedirs(dname)
def tointex(f):
    return int(f*100)
for i,plate in enumerate(plates):
    if len(plate) < 2: continue
    p0 = plate[0]
    # [print(f'{p.x-p0.x:5.2f}:{p.y-p0.y:5.2f} | ', end ='') for p in plate]
    # print()
    # print('-'*80)
    filenm=os.path.join(dname, f'{p0.name}.svg'.replace('>','v'))
    prevP = plate[-1]
    elements = []
    maxX = max([p.x for p in plate])
    maxY = max([p.y for p in plate])
    minX = min([p.x for p in plate])
    minY = min([p.y for p in plate])
    for p in plate:
        elements.append(
            svg.Line(
                x1=tointex(prevP.x-minX), x2=tointex(p.x-minX),
                y1=tointex(prevP.y-minY), y2=tointex(p.y-minY),
                stroke="blue",
                stroke_width=5,
            )
        )
        prevP = p
    canvas = svg.SVG(
        width=tointex((maxX - minX)*1.02),
        height=tointex((maxY - minY)*1.02),
        elements=elements,
    )
    with open(filenm, 'w') as f:
        f.write(canvas.__str__())
print('='*80)

# for br in branches:
#     points = br.getPoints()
#     points.append(points[0])
#     xs = np.array([p.x for p in points])
#     ys = np.array([p.y for p in points])
#     plt.plot(xs, ys)
#     plt.text((br.node2.point.x+br.node1.point.x)/2 , (br.node2.point.y+br.node1.point.y)/2, br.name)
#     # plate_ps = br.platePoints
#     # # plate_ps.append([0])
#     # plt.plot([p.x for p in plate_ps], [p.y for p in plate_ps], 'b-')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()



# V check extreme angles
# X produce production plot for each branch. length, 2angles, mark extreeme angle.
# X save general plot, and foreach branch plots to disk
# V list of branches in text: name, len, 2angs(mark extreme with "!")

# V calc base plates
# to svg base plates