import numpy as np
import networkx as nx

# Parsing from Desmos
radii_string = r"\left[1,0.5,0.5,0.5,0.5,0.3,1,0.8\right]"
centers_string = r"\left[\left(0,0\right),\left(0,1\right),\left(0,-1\right),\left(1,0\right),\left(-1,0\right),\left(0.75,0.5\right),\left(2.5,0\right),\left(2,1\right)\right]"
r = eval(radii_string.replace(r"\left[", "[").replace(r"\right]", "]"))
c = [np.array(i) for i in eval(centers_string.replace(r"\left[", "[").replace(r"\right]", "]").replace(r"\left(", "(").replace(r"\right)", ")"))]
circles = [[c[i], r[i]] for i in range(len(r))]

# circles = [[np.array(x,y), r],...]

NUM_OF_CIRCLES = len(circles)


def find_intersections(circ1: list[np.array, float], circ2: list[np.array, float]) -> list[np.array]:
    pos1 = circ1[0]
    pos2 = circ2[0]
    r1 = circ1[1]
    r2 = circ2[1]
    dist = np.linalg.norm(pos2 - pos1)
    if not (r1 - r2) ** 2 <= dist ** 2 <= (r1 + r2) ** 2:
        return [None, None]
    v = (pos2 - pos1) / dist
    a = 0.5 * (dist ** 2 + r1 ** 2 - r2 ** 2) / dist
    h = np.sqrt(r1 ** 2 - a ** 2)
    pt1 = pos1 + v * a + np.array([v[1], -v[0]]) * h
    pt2 = pos1 + v * a - np.array([v[1], -v[0]]) * h
    if np.array_equal(pt1, pt2):
        return [None, None]
    return [pt1, pt2]


def validate_point(circs, p, c1, c2=None):
    for index, circ in enumerate(circs):
        if np.linalg.norm(p - circ[0]) < circ[1] and index not in [c1, c2]:
            return True
    return False


intersections = {i: [] for i in range(NUM_OF_CIRCLES)}
int_points = []
pops = 0
for i in range(NUM_OF_CIRCLES):
    for j in range(i + 1, NUM_OF_CIRCLES):
        pts = list({array.tobytes(): array for array in
                    list(filter(lambda item: item is not None, find_intersections(circles[i], circles[j])))}.values())
        dummy = pts.copy()
        for ind, pt in enumerate(dummy):
            if validate_point(circles, pt, i, j):
                pts.pop(ind - pops)
                pops += 1
        pops = 0
        int_points += pts
        intersections[i] += pts
        intersections[j] += pts

node_value_dict = {ind: i for ind, i in enumerate([i[0] for i in circles] + int_points)}
value_node_dict = {v.tobytes(): k for k, v in node_value_dict.items()}

network = {k: [value_node_dict[i.tobytes()] for i in v] for k, v in intersections.items()}

for k, vlist in list(network.items()):
    for v in vlist:
        if v in network.keys():
            network[v].append(k)
        else:
            network[v] = [k]

G = nx.Graph(network)
isolated_circles = list(nx.isolates(G))
outer_isolated_circles = [k for k in isolated_circles if not validate_point(circles, circles[k][0], k)]
G.remove_nodes_from(isolated_circles)


def arc_angle(v2):
    return np.mod(np.arctan2(v2[1], v2[0]) + 2*np.pi, 2*np.pi)


total_sector_area = 0
for ind, i in enumerate(circles):
    connections = [(node_value_dict[k], k) for k in network[ind]]
    center = i[0]
    radius = i[1]
    angles = [(arc_angle(point[0] - center), point[0], point[1]) for point in connections]
    sorted_angles = sorted(angles)
    segments = []
    for j in range(len(sorted_angles)):
        segments.append((sorted_angles[j], sorted_angles[(j+1) % len(sorted_angles)]))

    arc_midpoints = [(arc[0][0] + arc[1][0])/2 for arc in segments]
    if len(arc_midpoints) == 2:
        arc_midpoints[1] = arc_midpoints[1] + np.pi
    exterior_arcs = [not validate_point(circles, np.array([np.cos(t), np.sin(t)])*radius + center, ind) for t in arc_midpoints]
    for n, arc in enumerate(segments):
        t1 = arc[0][0]
        t2 = arc[1][0]
        if exterior_arcs[n]:
            total_sector_area += 0.5 * radius**2 * ((t2-t1+2*np.pi) % (2*np.pi))
for node in outer_isolated_circles:
    total_sector_area += np.pi * circles[node][1]**2


def generate_sub_polys(graph: nx.Graph):
    sub_polys = []
    poly_crosses = [i[0] for i in graph.degree if i[1] > 2]
    if not poly_crosses:
        p = list(graph.nodes)[0]
        sub_poly = [p]
        dynpoint = list(graph.neighbors(p))[0]
        old_dynpoint = p
        while dynpoint != p:
            sub_poly.append(dynpoint)
            path = nx.to_dict_of_lists(graph)[dynpoint]
            path.remove(old_dynpoint)
            old_dynpoint = dynpoint
            dynpoint = path[0]
        return [sub_poly]
    for p in poly_crosses:
        nbors = graph.neighbors(p)
        for nbor in nbors:
            sub_poly = [p]
            dynpoint = nbor
            old_dynpoint = p
            while dynpoint != p:
                sub_poly.append(dynpoint)
                path = nx.to_dict_of_lists(graph)[dynpoint]
                path.remove(old_dynpoint)
                old_dynpoint = dynpoint
                dynpoint = path[0]
            if set(sub_poly) not in [set(i) for i in sub_polys]:
                sub_polys.append(sub_poly)
    return sub_polys


def get_subgraphs(graph: nx.Graph):
    return [graph.subgraph(comp).copy() for comp in nx.connected_components(graph)]


polygons = []
for g in get_subgraphs(G):
    polygons += [[node_value_dict[j] for j in i] for i in generate_sub_polys(g)]


def area(p):
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(p)))


def segments(p):
    return zip(p, p[1:] + [p[0]])


total_poly_area = sum([area(poly) for poly in polygons])
total_area = total_poly_area + total_sector_area
print("Total area:", total_area)
