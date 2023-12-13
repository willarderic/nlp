from os import link
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
import math
from re import findall
from sympy import *

RESIDUE_NAMES = ["ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "CYS", "SEC", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

protein_name = "3nir"

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
        print(list(res.get_atoms())[0])
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
                        orientation = '+'
                        if getLineAngle(line1) > math.pi:
                            orientation = '-'
                        bond_to_id[bond] = counter
                        gauss_code += 'B' + str(counter) + orientation
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
                        line_start_angle = getLineAngle(bond_start_line)
                        dot_prod = x1 * x2  + y1 * y2
                        if dot_prod >= 0 and line_start_angle <= math.pi:
                            gauss_code += "B" + str(counter) + "+"
                        elif dot_prod >= 0 and line_start_angle > math.pi:
                            gauss_code += "B" + str(counter) + "-"
                        elif dot_prod < 0 and line_start_angle <= math.pi:
                            gauss_code += "B" + str(counter) + "-"
                        else:
                            gauss_code += "B" + str(counter) + "+"
                        counter += 1
                    elif bond in bond_to_id.keys() and bond_to_id[bond] != None:
                        bond_start_line = lines[bond[0]]
                        x1, y1 = line1[1][0] - line1[0][0], line1[1][1] - line1[0][1]
                        x2, y2 = bond_start_line[1][0] - bond_start_line[0][0], bond_start_line[1][1] - bond_start_line[0][1]
                        line_start_angle = getLineAngle(bond_start_line)
                        dot_prod = x1 * x2  + y1 * y2
                        if dot_prod >= 0 and line_start_angle <= math.pi:
                            gauss_code += "B" + str(bond_to_id[bond]) + "+"
                        elif dot_prod >= 0 and line_start_angle > math.pi:
                            gauss_code += "B" + str(bond_to_id[bond]) + "-"
                        elif dot_prod < 0 and line_start_angle <= math.pi:
                            gauss_code += "B" + str(bond_to_id[bond]) + "-"
                        else:
                            gauss_code += "B" + str(bond_to_id[bond]) + "+"
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
    add_equals_y = False
    symbol_table = {}
    gauss_code = gauss_code[1:len(gauss_code)-1]

    # When index gets above 9, cutting every 3 characters doesnt work  
    #arc_ends = [(gauss_code[i:i+3]) for i in range(0, len(gauss_code), 3)]
    # Use regex instead
    arc_ends = findall(r'[B|O|U]{1}\d+[-|+]{1}', gauss_code)
    # Create a map of arc ends so that I can find whether a bond is positive or negative
    arc_end_map = {}
    bond_tracker = {}
    
    x_start = None
    y_start = None
    for i, arc_end in enumerate(arc_ends):
        arc_idx = arc_end[1]
        if arc_idx == '1' and x_start == None:
            x_start = i
        elif arc_idx == '1':
            y_start = i
        
        if arc_end[0] == 'B' and not bond_tracker.get(arc_end):
            bond_tracker[arc_end[:2]] = False

        found_key = None
        for key in arc_end_map:
            if arc_idx == key[0][1]:
                found_key = key
        
        if found_key == None:
            arc_end_map[(arc_end, i)] = None
        else:
            arc_end_map[(arc_end, i)] = found_key
            arc_end_map[found_key] = (arc_end, i)

    # print(arc_end_map)
    # print(x_start, y_start)
    n = 0
    expr = 'x'
    swp_expr = 'y'
    idx = x_start
    swp_idx = y_start
    while n < len(arc_ends) + 1:
        # print(n, idx, swp_idx, symbol_table)
        if idx >= len(arc_ends) or idx == swp_idx:
            tmp = idx
            idx = swp_idx
            swp_idx = tmp

            tmp = expr
            expr = swp_expr
            swp_expr = tmp

        arc_end = arc_ends[idx]
        # Check if we need to swap
        if arc_end[0] == 'U':
            if symbol_table.get(arc_end_map[(arc_end, idx)]) == None:
                tmp = idx
                idx = swp_idx
                swp_idx = tmp

                tmp = expr
                expr = swp_expr
                swp_expr = tmp
        elif arc_end[0] == 'B':
            _, other_idx = arc_end_map[(arc_end, idx)]
            if (other_idx != x_start and other_idx != y_start) and symbol_table.get((arc_ends[other_idx - 1], other_idx - 1)) == None:
                tmp = idx
                idx = swp_idx
                swp_idx = tmp

                tmp = expr
                expr = swp_expr
                swp_expr = tmp
        arc_end = arc_ends[idx]

        if arc_end[0] == 'U':
            if arc_end[2] == '+':
                op = 'op'
            elif arc_end[2] == '-':
                op = 'inv'

            expr = f"{expr} {symbol_table[arc_end_map[(arc_end, idx)]]} {op}"
        elif arc_end[0] == 'B':
            linked_arc_end, linked_arc_end_idx = arc_end_map[(arc_end, idx)]
            swap_params = False
            if linked_arc_end_idx > idx:
                if arc_end[2] == '+' and linked_arc_end[2] == '+':
                    op = 'R2'
                    swap_params = True
                elif arc_end[2] == '-' and linked_arc_end[2] == '-':
                    op = 'R1'
                    swap_params = False
                else:
                    op = 'R3'
                    swap_params = False
                _, other_idx = arc_end_map[(arc_end, idx)]
                operand = ''
                if other_idx == y_start:
                    operand = 'y'
                elif other_idx == x_start:
                    operand = 'x'
                else:
                    operand = symbol_table[(arc_ends[other_idx - 1], other_idx - 1)]
                if swap_params:
                    expr = f"{operand} {expr} {op}"
                else:
                    expr = f"{expr} {operand} {op}"
            elif linked_arc_end_idx < idx:
                if arc_end[2] == '+' and linked_arc_end[2] == '+':
                    op = 'R1'
                    swap_params = False
                elif arc_end[2] == '-' and linked_arc_end[2] == '-':
                    op = 'R2'
                    swap_params = True
                else:
                    op = 'R3'
                    swap_params = False
                _, other_idx = arc_end_map[(arc_end, idx)]
                operand = ''
                if other_idx == y_start:
                    operand = 'y'
                elif other_idx == x_start:
                    operand = 'x'
                else:
                    operand = symbol_table[(arc_ends[other_idx - 1], other_idx - 1)]
                if swap_params:
                    expr = f"{operand} {expr} {op}"
                else:
                    expr = f"{expr} {operand} {op}"
        
        # Handles the over crossing as well
        symbol_table[(arc_end, idx)] = expr
        idx += 1
        n += 1
    return symbol_table[arc_ends[y_start - 1], y_start - 1]
print(f"y={gen_algebra('NB1+U2-B1+O2-C')}")
print(f"y={gen_algebra('NB1+O2-B1+U2-C')}")
print(f"y={gen_algebra('NB1+O2-B3+U4-O5-U6-B1-O6-U2-B3+O4-U5-C')}")
print(f"y={gen_algebra('NB1+U2-B3+O4-U5-O6-B1-U6-O2-B3+U4-O5-C')}")

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        return None
    else:
        return x % m

def prerocess_var(s, x, y):
    if s == 'x':
        return x
    if s == 'y':
        return y
    return s

# function to evaluate reverse polish notation 
def evaluate(expression, a_in, b_in, n, m_in):
    print(expression)
    # A = a^-1
    x, y, a, b, m, A = symbols("x y a b m A")

    op = a * x + (1 - a) * y
    inv =  A * x + (1 - A) * y
    R1 = b * x + (1 - b) * y
    R2 = a * (1 - b) * x + (b + (1 - a) * (1 - b)) * y
    R3 = m * x + (1 - m) * y

    op = op.subs(a, a_in)
    inv = inv.subs(A, modinv(a_in, n))
    R1 = R1.subs(b, b_in)
    R2 = R2.subs([(a, a_in), (b, b_in)], simultaneous=True)
    print(R2)
    R3 = R3.subs(m, m_in)

    expression = expression.split() 
    expression = [prerocess_var(e, x, y) for e in expression]
    stack = [] 
        
    # RPN Evaluation
    for ele in expression:
        if ele not in ['op', 'inv', 'R1', 'R2', 'R3']:
            stack.append(ele)
        else:
            right = stack.pop() 
            left = stack.pop()

            if ele == 'op':
                stack.append(op.subs([(x, left), (y, right)], simultaneous=True))
            elif ele == 'inv':
                stack.append(inv.subs([(x, left), (y, right)], simultaneous=True))
            elif ele == 'R1': 
                stack.append(R1.subs([(x, left), (y, right)], simultaneous=True)) 
            elif ele == 'R2': 
                stack.append(R2.subs([(x, left), (y, right)], simultaneous=True)) 
            elif ele == 'R3':
                stack.append(R3.subs([(x, left), (y, right)], simultaneous=True))
        
    # return final answer. 
    return stack.pop() 

def calc_colorings(expression, a, b, n, m):
    x, y = symbols('x y')
    equation = evaluate(expression, a, b, n, m)
    print(f'Equation={equation}')
    count = 0
    for i in range(n):
        for j in range(n):
            result = equation.subs([(x, i), (y, j)])
            equals = True if (result % n) == j else False
            # print(f'{j}={result % n} is {equals}')
            if equals:
                count += 1

    print('Number of colorings =', count)

calc_colorings(gen_algebra('NB1+U2-B1+O2-C'), 8, 2, 15, 6)
calc_colorings(gen_algebra('NB1+O2-B1+U2-C'), 8, 2, 15, 6)
calc_colorings(gen_algebra('NB1+O2-B3+U4-O5-U6-B1-O6-U2-B3+O4-U5-C'), 7, 8, 15, 6)
calc_colorings(gen_algebra(gauss_code), 8, 2, 15, 6)
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
# plt.text(ys[0], zs[0], "N")
# plt.text(ys[-1], zs[-1], "C")
# for bl in bond_lines:
#     plt.plot([bl[0][0], bl[1][0]], [bl[0][1], bl[1][1]], color="purple")

# plt.show()



