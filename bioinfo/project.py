from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
import math

RESIDUE_NAMES = ["ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "CYS", "SEC", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

protein_name = "5awl"

structure_id = protein_name
filename = f"/Users/ericwillard/School/bioinfo/project/pdb/{protein_name}.pdb"
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(structure_id, filename)[0]
xs = []
ys = []
zs = []
residues = []

def parse_ss_bonds(filename):
    with open(filename) as f:
        return [(int(line[17:21].strip()), int(line[31:35].strip())) for line in f.readlines() if line[:6] == "SSBOND"]

ss_bonds = parse_ss_bonds(filename)
print(ss_bonds)
print(f"Number of residues = {len(list(structure.get_residues()))}")
for res in structure.get_residues():
    # Change this from blacklist to whitelist (only accept the 20 amino acids, nothing else)
    if res.get_resname() in RESIDUE_NAMES:
        residues.append(res.get_resname())

        x, y, z = list(res.get_atoms())[0].get_coord()
        xs.append(x)
        ys.append(y)
        zs.append(z)

class Intersection:
    counter = 1

    def __init__(self, x, y):
        self.y = y
        self.z = z
        self.id = Intersection.counter
        Intersection.counter += 1

class Segment:
    counter = 1

    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.id = Segment.counter
        Segment.counter += 1

    def intersect2d(self, segment):
        x1 = self.p1[1]
        y1 = self.p1[2]
        x2 = self.p2[1]
        y2 = self.p2[2]
        x3 = segment.p1[1]
        y3 = segment.p1[2]
        x4 = segment.p2[1]
        y4 = segment.p2[2]

        if (x2 - x1 == 0) or (x4 - x3 == 0):
            return False
        
        m1 = (y2 - y1) / (x2 - x1)
        m2 = (y4 - y3) / (x4 - x3)

        # same line or parallel
        if (m1 == m2):
            return False

        b1 = y1 - m1 * x1
        b2 = y3 - m2 * x3

        if (m1 - m2) == 0:
            return False

        x_int = (b2 - b1) / (m1 - m2)

        if (x_int < max( min(x1,x2), min(x3,x4))) or (x_int > min( max(x1,x2), max(x3,x4))):
            return False
        else:
            inter = Intersection(x_int, b1 + m1 * x_int)
            return inter
        
    def calculate_x(self, intersection2d):
        x1, y1, z1 = self.p1[0], self.p1[1], self.p1[2]
        x2, y2, z2 = self.p2[0], self.p2[1], self.p2[2]
        print(x1, y1, z1)
        x = ((x1 - x2) * (intersection2d.z - z1) / (z1 - z2)) + x1
        return x

    def getLineAngle(self):
        x, y = self.p2[0] - self.p1[0], self.p2[1] - self.p1[1]
        if x == 0:
            if y > 0:
                return math.pi
            elif y < 0:
                return 3 * math.pi / 2

        theta = math.atan(y/x)
        return theta


points3d = list(zip(xs, ys, zs))
segments = []
for i in range(len(points3d) - 1):
    segments.append(Segment(points3d[i], points3d[i+1]))
    print(f"{segments[i].id} -> ({segments[i].p1}, {segments[i].p2})")

intersections = []
gauss_code = 'N'
for i, line1 in enumerate(segments):
    for j, line2 in enumerate(segments): 
        if i != j and i != j - 1 and i != j + 1:
            if intersection := line1.intersect2d(line2):
                # Find over/under
                x1 = line1.calculate_x(intersection)
                x2 = line2.calculate_x(intersection)

                theta1 = line1.getLineAngle()
                theta2 = line2.getLineAngle()

                if x1 > x2:
                    gauss_code += 'O'
                    print("Line 1 over line 2")
                    if theta2 < theta1 and theta2 > (theta1 - math.pi):
                        gauss_code += '-'
                        print("negative crossing")
                    else:
                        gauss_code += '+'
                        print("positive crossing")
                else:
                    gauss_code += 'U'
                    print("Line 2 over line 1")
                    if theta2 < theta1 and theta2 > (theta1 - math.pi):
                        gauss_code += '+'
                        print("positive crossing")
                    else:
                        gauss_code += '-'
                        print("negative crossing")
                
                # print(f"Line {lines3d[i]} with intersection {intersection} and x = {x1}")
                # print(f"Line {lines3d[j]} with intersection {intersection} and x = {x2}")
                # fig =plt.figure()
                # ax = fig.add_subplot(111, projection='3d')
                # ax.plot([lines3d[i][0][0] ,lines3d[i][1][0]], [lines3d[i][0][1] ,lines3d[i][1][1]], [lines3d[i][0][2] ,lines3d[i][1][2]], color='blue')
                # ax.plot([lines3d[j][0][0] ,lines3d[j][1][0]], [lines3d[j][0][1] ,lines3d[j][1][1]], [lines3d[j][0][2] ,lines3d[j][1][2]], color='red')
                # ax.scatter(x1, intersection[0], intersection[1], color='green', s=50)
                # ax.scatter(x2, intersection[0], intersection[1], color='green', s=50)
                # plt.show()

                # Find 
                
                intersections.append((intersection.y, intersection.z))

gauss_code += 'C'
print(gauss_code)

plt.plot(ys, zs)
plt.scatter(*zip(*intersections), color='red')

plt.show()
