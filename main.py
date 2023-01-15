import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib import patches as ptc


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


def generate_network(circs):
    intersections = {i: [] for i in range(NUM_OF_CIRCLES)}
    int_points = []
    pops = 0
    for i in range(NUM_OF_CIRCLES):
        for j in range(i + 1, NUM_OF_CIRCLES):
            pts = list({array.tobytes(): array for array in
                        list(filter(lambda item: item is not None, find_intersections(circs[i], circs[j])))}.values())
            dummy = pts.copy()
            for ind, pt in enumerate(dummy):
                if validate_point(circs, pt, i, j):
                    pts.pop(ind - pops)
                    pops += 1
            pops = 0
            int_points += pts
            intersections[i] += pts
            intersections[j] += pts

    node_value_dict = {ind: i for ind, i in enumerate([i[0] for i in circs] + int_points)}
    value_node_dict = {v.tobytes(): k for k, v in node_value_dict.items()}

    network = {k: [value_node_dict[i.tobytes()] for i in v] for k, v in intersections.items()}

    for k, vlist in list(network.items()):
        for v in vlist:
            if v in network.keys():
                network[v].append(k)
            else:
                network[v] = [k]

    graph = nx.Graph(network)
    isolated_circles = list(nx.isolates(graph))
    outer_isolated_circles = [k for k in isolated_circles if not validate_point(circs, circs[k][0], k)]

    return graph, node_value_dict, outer_isolated_circles


def arc_angle(v2):
    return np.mod(np.arctan2(v2[1], v2[0]) + 2 * np.pi, 2 * np.pi)


def area_of_sectors(circs, graph, node_value_dict, outer_isolated_circles):
    total_sector_area = 0
    outer_arcs = {}
    for ind, i in enumerate(circs):
        connections = [(node_value_dict[k], k) for k in nx.neighbors(graph, ind)]
        center = i[0]
        radius = i[1]
        angles = [(arc_angle(point[0] - center), point[0], point[1]) for point in connections]
        sorted_angles = sorted(angles)
        segments = []
        for j in range(len(sorted_angles)):
            segments.append((sorted_angles[j], sorted_angles[(j + 1) % len(sorted_angles)]))

        arc_midpoints = [(arc[0][0] + arc[1][0]) / 2 if arc[0][0] < arc[1][0] else (arc[0][0] + arc[1][0]) / 2 + np.pi
                         for arc in segments]
        exterior_arcs = [not validate_point(circs, np.array([np.cos(t), np.sin(t)]) * radius + center, ind) for t in
                         arc_midpoints]
        for n, arc in enumerate(segments):
            t1 = arc[0][0]
            t2 = arc[1][0]
            if exterior_arcs[n]:
                outer_arcs[arc[0][2]] = [ind, arc[1][2]]
                total_sector_area += 0.5 * radius ** 2 * ((t2 - t1 + 2 * np.pi) % (2 * np.pi))
    for node in outer_isolated_circles:
        total_sector_area += np.pi * circs[node][1] ** 2
    return total_sector_area, outer_arcs


def get_subgraphs(graph: nx.Graph):
    return [graph.subgraph(comp).copy() for comp in nx.connected_components(graph)]


def area_of_polygons(maingraph: nx.Graph, node_value_dict: dict, outer_arcs: dict):
    subgraphs: list[nx.Graph] = get_subgraphs(maingraph)
    total_area_polygons = 0
    for graph in subgraphs:
        if graph.number_of_nodes() == 1:
            continue
        areas = []
        couples = {}
        graphnodes = list(graph.nodes)
        for start, couple in outer_arcs.items():
            if couple[1] in graphnodes:
                couples[start] = couple
        visited = []
        startnode = None
        for n in graph.nodes:
            if n > NUM_OF_CIRCLES:
                startnode = n
                break
        while not sorted(visited) == sorted(graphnodes):
            path = [startnode]
            while path[-1] != startnode or len(path) == 1:
                c_from_dn = couples[path[-1]]
                for n in c_from_dn:
                    path.append(n)
            path.pop(-1)
            points = []
            for n in path:
                points.append(node_value_dict[n])
                if n not in visited:
                    visited.append(n)
            areas.append(area(points))
            for n in graphnodes:
                if n not in visited and n > NUM_OF_CIRCLES:
                    startnode = n
                    break
        area_of_subgraph = max(areas)
        areas.remove(area_of_subgraph)
        total_area_polygons += area_of_subgraph - sum(areas)
    return total_area_polygons


def area(p):
    def segments(pt):
        return zip(pt, pt[1:] + [pt[0]])

    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for ((x0, y0), (x1, y1)) in segments(p)))


def get_total_area(circs):
    graph, node_value_dict, outer_isolated_circles = generate_network(circs)
    a_s, arcs = area_of_sectors(circs, graph, node_value_dict, outer_isolated_circles)
    a_p = area_of_polygons(graph, node_value_dict, arcs)
    return a_s + a_p, graph, node_value_dict


def draw_graph(circs, graph, node_value_dict):
    fig, ax = plt.subplots()
    for c in circs:
        ax.add_patch(ptc.Circle(xy=c[0], radius=c[1],
                                fill=True, linewidth=None,
                                fc=(0, 0.4, 0, 0.2)))
    pos = {k: node_value_dict[k] for k in list(g.nodes)}
    nx.draw(graph, pos=pos,
            node_size=[50]*NUM_OF_CIRCLES+[10]*(graph.number_of_nodes()-NUM_OF_CIRCLES),
            node_color=["red"]*NUM_OF_CIRCLES+["blue"]*(graph.number_of_nodes()-NUM_OF_CIRCLES))
    nx.draw_networkx_labels(graph, pos=pos, font_size=6,
                            labels={i: str(i) for i in range(NUM_OF_CIRCLES)})
    ax.set_aspect('equal', adjustable='box')
    plt.show()


circles = [
    [np.array([3560, -92]), 35],
    [np.array([3527, -80]), 37],
    [np.array([3495, -102]), 40],
    [np.array([3523, -132]), 43],
    [np.array([3563, -149]), 46],
    [np.array([3620, -209]), 50],
    [np.array([3682, -260]), 54],
    [np.array([3650, -92]), 30],
    [np.array([3610, -65]), 30],
    [np.array([3620, -125]), 20],
    [np.array([3400, -50]), 20],
    [np.array([3450, -40]), 20],
    [np.array([3430, -70]), 20],
    [np.array([3420, -20]), 20],
    [np.array([3500, -240]), 40]
]

NUM_OF_CIRCLES = len(circles)

A, g, nvd = get_total_area(circles)
draw_graph(circles, g, nvd)
print("Total area:", A)
