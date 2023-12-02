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

# # Draw 3d shape
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.plot(xs, ys, zs, color='red')
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')

# for x, y, z, residue in zip(xs, ys, zs, residues):
#     ax.text(x, y, z, residue)

# plt.show()

def line_intersection(line1, line2):
    x1 = line1[0][0]
    y1 = line1[0][1]
    x2 = line1[1][0]
    y2 = line1[1][1]
    x3 = line2[0][0]
    y3 = line2[0][1]
    x4 = line2[1][0]
    y4 = line2[1][1]

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
        return (x_int, b1 + m1 * x_int)


points = list(zip(ys, zs))
# print("Points", points)
lines = []
for i in range(len(points) - 1):
    lines.append((points[i], points[i+1]))

points3d = list(zip(xs, ys, zs))
lines3d = []
for i in range(len(points3d) - 1):
    lines3d.append((points3d[i], points3d[i+1]))
# print(lines3d)

def calculate_x(line3d, intersection2d):
    x1, y1, z1 = line3d[0]
    x2, y2, z2 = line3d[1]
    # print(x1, y1, z1)
    x = ((x1 - x2) * (intersection2d[1] - z1) / (z1 - z2)) + x1
    return x

def getLineAngle(line):
    x, y = line[1][0] - line[0][0], line[1][1] - line[0][1]
    # print(x, y)
    if x >  0 and y > 0:
        # QUAD 1
        theta = math.atan(abs(y)/abs(x))
        return theta
    elif x < 0 and y > 0:
        # QUAD 2
        theta = math.atan(abs(y)/abs(x))
        # atan returns the negative angle from straight up here, so if we make it positive and add 90 degrees, we get the angle from (0,1)
        return theta + math.pi / 2
    elif x < 0 and y < 0:
        # QUAD 3
        theta = math.atan(abs(y)/abs(x))
        return theta + math.pi
    elif x > 0 and y < 0:
        # QUAD 4
        theta = math.atan(abs(y)/abs(x))
        return theta + 3 * math.pi / 2
    elif x == 0 and y > 0:
        return math.pi / 2
    elif x == 0 and y < 0:
        return 3 * math.pi / 2
    elif x > 0 and y == 0:
        return 0
    elif x < 0 and y == 0:
        return math.pi

# print("Line", lines)
intersections = []

gauss_code = 'N'
intersection_to_id = {}
bond_to_id = {}
counter = 1
for i, line1 in enumerate(lines):
    for j, line2 in enumerate(lines): 
        if i != j and i != j - 1 and i != j + 1:
            # currently not handling the orientation of the bond. 
            # Can change the '+' on the first B gauss code to a '-' if y < 0 and keep a '+' if y > 0
            bond = None
            if i in [ss_bonds[k][0] for k in range(len(ss_bonds))]:
                for k in range(len(ss_bonds)):
                    if i == ss_bonds[k][0]:
                        bond = ss_bonds[k]
                if bond != None:
                    if not bond in bond_to_id.keys():
                        bond_to_id[bond] = counter
                        gauss_code += 'B' + str(counter) + '+'
                        counter += 1

            bond = None
            if i in [ss_bonds[k][1] for k in range(len(ss_bonds))]:
                for k in range(len(ss_bonds)):
                    if i == ss_bonds[k][1]:
                        bond = ss_bonds[k]
                if bond != None:
                    if not bond in bond_to_id.keys():
                        bond_to_id[bond] = counter
                        bond_start_line = lines[bond[0]]
                        x1, y1 = line1[1][0] - line1[0][0], line1[1][1] - line1[0][1]
                        x2, y2 = bond_start_line[1][0] - bond_start_line[0][0], bond_start_line[1][1] - bond_start_line[0][1]
                        
                        dot_prod = x1 * x2  + y1 * y2
                        print(f"<({x1}, {y1}), ({x2}, {y2})> = {dot_prod}")
                        if dot_prod >= 0:
                            gauss_code += "B" + str(counter) + "+"
                        else:
                            gauss_code += "B" + str(counter) + "-"
                        counter += 1
                    elif bond in bond_to_id.keys() and bond_to_id[bond] != None:
                        bond_start_line = lines[bond[0]]
                        x1, y1 = line1[1][0] - line1[0][0], line1[1][1] - line1[0][1]
                        x2, y2 = bond_start_line[1][0] - bond_start_line[0][0], bond_start_line[1][1] - bond_start_line[0][1]
                        
                        dot_prod = x1 * x2  + y1 * y2
                        print(f"<({x1}, {y1}), ({x2}, {y2})> = {dot_prod}")
                        if dot_prod >= 0:
                            gauss_code += "B" + str(bond_to_id[bond]) + "+"
                        else:
                            gauss_code += "B" + str(bond_to_id[bond]) + "-"
                        bond_to_id[bond] = None

            if intersection := line_intersection(line1, line2):
                intersections.append(intersection)
                if (not (i, j) in intersection_to_id.keys()) and (not (j, i) in intersection_to_id.keys()):
                    intersection_to_id[(i, j)] = counter
                    counter += 1
                
                # Find over/under
                x1 = calculate_x(lines3d[i], intersection)
                x2 = calculate_x(lines3d[j], intersection)

                theta1 = getLineAngle(line1)
                theta2 = getLineAngle(line2)
                # print(f"Theta1 = {theta1}, theta2 = {theta2}")

                inter_num = intersection_to_id[(i, j)] if (i, j) in intersection_to_id.keys() else intersection_to_id[(j, i)]

                if x1 > x2:
                    gauss_code += 'O' + str(inter_num)
                    # print("Line 1 over line 2")
                    if theta2 < theta1 and theta2 > (theta1 - math.pi):
                        gauss_code += '-'
                        # print("negative crossing")
                    else:
                        gauss_code += '+'
                        # print("positive crossing")
                else:
                    gauss_code += 'U' + str(inter_num)
                    # print("Line 2 over line 1")
                    if theta2 < theta1 and theta2 > (theta1 - math.pi):
                        gauss_code += '+'
                        # print("positive crossing")
                    else:
                        gauss_code += '-'
                        # print("negative crossing")
                
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
                
                
# print(intersection_to_id)
gauss_code += 'C'
print(gauss_code)

def gen_algebra(gauss_code):
    symbol_table = {}
    gauss_code = gauss_code[1:len(gauss_code)-1]
    arc_ends = [(gauss_code[i:i+3]) for i in range(0, len(gauss_code), 3)]
    # Create a map of arc ends so that I can find whether a bond is positive or negative
    arc_end_map = {}
    for arc_end in arc_ends:
        idx = arc_end[1]

        found_key = None
        for key in arc_end_map:
            if idx == key[1]:
                found_key = key
        
        if found_key == None:
            arc_end_map[arc_end] = None
        else:
            arc_end_map[arc_end] = found_key
            arc_end_map[found_key] = arc_end

    expr = 'x'
    expr_list = []

    def next_var():
        next_var.var_cnt += 1
        return f'x{next_var.var_cnt}'
    next_var.var_cnt = 0

    for arc_end in arc_ends:
        # Save the resulting expression from the arc
        n_var = next_var()
        expr_list.append(f'{n_var}={expr}')
        symbol_table[arc_end] = n_var

        if arc_end[0] == 'U':
            op = ''
            if arc_end[2] == '+':
                op = 'op'
            elif arc_end[2] == '-':
                op = 'inv'
            
            if var := symbol_table.get(arc_end_map[arc_end]):
                expr = f'{op}({expr},{var})'
            else:
                n_var = next_var()
                symbol_table[arc_end_map[arc_end]] = n_var
                expr = f'{op}({expr},{n_var})'
                
        elif arc_end[0] == 'B':

            op1 = ''
            op2 = ''
            linked_arc_end = arc_end_map[arc_end]
            if arc_end[2] == '+' and linked_arc_end[2] == '+':
                op1 = 'R1'
                op2 = 'R2'
            elif arc_end[2] == '-' and linked_arc_end[2] == '-':
                op1 = 'R2'
                op2 = 'R1'
            else:
                op1 = 'R3'
                op2 = 'R3'
            
            if var := symbol_table.get(arc_end_map[arc_end]):
                expr = f'{op}({expr},{var})'



        


gen_algebra(gauss_code)

# print(intersections)
# Draw 2D projection
# plt.clf()
# bond_lines = [(points[i - 1], points[j - 1]) for i, j in ss_bonds]
# print(bond_lines)

# fig =plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(xs, ys, zs, color='blue')
# # plt.scatter(*zip(*intersections), color='red')
# plt.show()


# plt.plot(ys, zs)
# plt.scatter(*zip(*intersections), color='red')
# for bl in bond_lines:
#     plt.plot([bl[0][0], bl[1][0]], [bl[0][1], bl[1][1]], color="purple")

# plt.show()



