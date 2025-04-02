#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last update December 20, 2023
"""
import numpy as np
import scipy.io as sio
import argparse
import os
#import matlab.engine %no need
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
import random
import heapq
import itertools
import time

def get_complete_graph_rod(n_grid=6):
    """
    Get the complete graph of a rod-based structure

    Rod-based structure:
        each vertex in the graph represents a node
        each edge in the graph represents a rod
        each node can be connnected to the 26 neighbors in the 3D space

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph
    # going through all the vertices (x,y,z)
    for x in range(n_grid):
        for y in range(n_grid):
            for z in range(n_grid):
                graph[(x,y,z)] = []
                # consider all neighbors, dx, dy, dz are the offsets between (x,y,z) and the neighbor (x_,y_,z_)
                for dx in [-1,0,1]:
                    for dy in [-1,0,1]:
                        for dz in [-1,0,1]:
                            if dx==dy==dz==0:
                                # skip self connection
                                continue
                            x_, y_, z_ = x+dx, y+dy, z+dz
                            if 0<=x_<n_grid and 0<=y_<n_grid and 0<=z_<n_grid:
                                # add the edges iff (x_,y_,z_) is not out of boundary
                                graph[(x,y,z)].append((x_,y_,z_))

    return graph

def get_complete_graph_plane(n_grid=6):
    """
    Get the complete graph of a plane-based structure, specifically triangle

    Plane-based structure:
        each vertex in the graph represents a center of a square
        each edge in the graph means two centers of squares are connected by two triangles connecting the centers and their shared edge
        each center of square are connected to the 12 neighboring square centers in the 3D space

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph

    # consider 3 cases because there are three orientations of the squares in the plane-based design

    # going through all the grids (x,y,z)
    for x in range(n_grid):
        for y in range(n_grid-1):
            for z in range(n_grid-1):
                x_c,y_c,z_c = x,y+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the square, the real vertices we consider
                graph[(x_c,y_c,z_c)] = []
                # consider all its neighbors
                for dx, dy, dz in [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(.5,.5,0),(-.5,.5,0),(.5,-.5,0),(-.5,-.5,0),(.5,0,.5),(-.5,0,.5),(.5,0,-.5),(-.5,0,-.5)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if 0<=x_<=n_grid-1 and 0<=y_<=n_grid-1 and 0<=z_<=n_grid-1:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))

    for y in range(n_grid):
        for x in range(n_grid-1):
            for z in range(n_grid-1):
                x_c,y_c,z_c = x+0.5,y,z+0.5
                graph[(x_c,y_c,z_c)] = []
                for dx, dy, dz in [(1,0,0),(-1,0,0),(0,0,1),(0,0,-1),(.5,.5,0),(-.5,.5,0),(.5,-.5,0),(-.5,-.5,0),(0,.5,.5),(0,-.5,.5),(0,.5,-.5),(0,-.5,-.5)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if 0<=x_<=n_grid-1 and 0<=y_<=n_grid-1 and 0<=z_<=n_grid-1:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))

    for z in range(n_grid):
        for x in range(n_grid-1):
            for y in range(n_grid-1):
                x_c,y_c,z_c = x+0.5,y+0.5,z
                graph[(x_c,y_c,z_c)] = []
                for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(.5,0,.5),(-.5,0,.5),(.5,0,-.5),(-.5,0,-.5),(0,.5,.5),(0,-.5,.5),(0,.5,-.5),(0,-.5,-.5)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if 0<=x_<=n_grid-1 and 0<=y_<=n_grid-1 and 0<=z_<=n_grid-1:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))


    return graph

def get_complete_graph_pyramid(n_grid=6):
    """
    Get the complete graph of a pyramid-based structure

    Pyramid-based structure:
        each vertex in the graph represents a center of a cube
        each edge in the graph means two centers of cubes are connected and share a squares
        each center of cubes are connected to the 6 neighboring cubes' centers in the 3D space

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph

    # going through all the grids (x,y,z)
    for x in range(-1,n_grid-1):
        for y in range(-1,n_grid-1):
            for z in range(-1,n_grid-1):
                #print("Node")
                #print(x,y,z)
                x_c,y_c,z_c = x+0.5,y+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
                # print("center node")
                # print(x_c,y_c,z_c)
                graph[(x_c,y_c,z_c)] = []
                # consider all its neighbors
                for dx, dy, dz in [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if -0.5<=x_<=n_grid-1 and -0.5<=y_<=n_grid-1 and -0.5<=z_<=n_grid-1:
                        # Check if the neighbor is within the boundary groups
                        if not (x_c in (-0.5, n_grid - 1.5) and x_ in (-0.5, n_grid - 1.5)):
                            if not (y_c in (-0.5, n_grid - 1.5) and y_ in (-0.5, n_grid - 1.5)):
                                if not (z_c in (-0.5, n_grid - 1.5) and z_ in (-0.5, n_grid - 1.5)):
                                    # Add the neighbor as a connection
                                    graph[(x_c, y_c, z_c)].append((x_, y_, z_))
                        #         else:
                        #             print(x_,y_,z_)
                        #     else:
                        #         print(x_,y_,z_)
                        # else:
                        #     print(x_,y_,z_)

                        # graph[(x_c,y_c,z_c)].append((x_, y_, z_))
                        # print("center node final")
                        # print(x_,y_,z_)
                    # else:
                    #     print(x_,y_,z_)


    return graph

def get_complete_graph_pyramid_old(n_grid=6):
    """
    Get the complete graph of a pyramid-based structure

    Plane-based structure:
        each vertex in the graph represents a center of a cube
        each edge in the graph means two centers of cubes are connected by two squares connecting the centers and their shared surface
        each center of cubes are connected to the 6 neighboring cubes' centers in the 3D space

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph

    # going through all the grids (x,y,z)
    for x in range(-1,n_grid-1):
        for y in range(-1,n_grid-1):
            for z in range(-1, n_grid-1):
                #print("Node")
                #print(x,y,z)
                x_c,y_c,z_c = x+0.5,y+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
                #print("center node")
                #print(x_c,y_c,z_c)
                graph[(x_c,y_c,z_c)] = []
                # consider all its neighbors
                for dx, dy, dz in [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if -0.5<=x_<=n_grid-1 and -0.5<=y_<=n_grid-1 and -0.5<=z_<=n_grid-1:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))
                        #print("center node final")
                        #print(x_,y_,z_)

    return graph

def get_complete_graph_cube(n_grid=6):
    """
    Get the complete graph of a cube-based structure

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph

    # going through all the grids (x,y,z)
    for x in range(n_grid):
        for y in range(n_grid):
            for z in range(n_grid):
                #print("Node")
                #print(x,y,z)
                x_c,y_c,z_c = x+0.5,y+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
                #print("center node")
                #print(x_c,y_c,z_c)
                graph[(x_c,y_c,z_c)] = []
                # consider all its neighbors
                for dx, dy, dz in [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if 0<x_<n_grid and 0<y_<n_grid and 0<z_<n_grid:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))
                        #print("center node final")
                        #print(x_,y_,z_)

    return graph

def get_interior_graph_cube(n_grid=6):
    """
    For periodic cube-base structure
    
    Get the interior graph of a cube-based structure

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph = {} # initialize the graph

    # going through all the grids (x,y,z)
    for x in range(1, n_grid-1):
        for y in range(1, n_grid-1):
            for z in range(1, n_grid-1):
                #print("Node")
                #print(x,y,z)
                x_c,y_c,z_c = x+0.5,y+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
                #print("center node")
                #print(x_c,y_c,z_c)
                graph[(x_c,y_c,z_c)] = []
                # consider all its neighbors
                for dx, dy, dz in [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)]:
                    x_, y_, z_ = x_c+dx, y_c+dy, z_c+dz
                    if 0.5<x_<n_grid-0.5 and 0.5<y_<n_grid-0.5 and 0.5<z_<n_grid-0.5:
                        graph[(x_c,y_c,z_c)].append((x_, y_, z_))
                        #print("center node final")
                        #print(x_,y_,z_)

    return graph

def get_exterior_graph_cube(n_grid=6):
    """
    For periodic cube-base structure
    
    Get the exterior graph of a cube-based structure

    input: n_grid, number of vertices on each side of the cube
    output: graph, the complete graph in the design space, represented as adjacency list
    """
    graph_xy_0 = {} # initialize the graph
    graph_yz_0 = {} # initialize the graph
    graph_xz_0 = {} # initialize the graph

    # going through all the grids for the boundaries at the origin

    z_0=0 #Surface xy at origine
    for x in range(1, n_grid-1):
        for y in range(1, n_grid-1):
            # print("Node")
            # print(x,y,z)
            x_c,y_c,z_0_c = x+0.5,y+0.5,z_0+0.5  # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
            # print("center node")
            # print(x_c,y_c,z_c)
            graph_xy_0[(x_c,y_c,z_0_c)] = []
            # consider all its neighbors
            for dx, dy in [(0,1),(0,-1),(1,0),(-1,0)]:
                x_, y_ = x_c+dx, y_c+dy
                if 0<x_<n_grid and 0<y_<n_grid:
                    graph_xy_0[(x_c,y_c,z_0_c)].append((x_, y_, z_0_c))
                    # print("center node final")
                    # print("XYZ0")
                    # print("XYZN")
                    # print(x_,y_,z_0_c)


    x_0=0 #Surface yz at origine
    #We do not redefine shared edges with xy
    for y in range(1, n_grid-1):
        for z in range(1, n_grid-1):
            # print("Node")
            # print(x,y,z)
            x_0_c,y_c,z_c = x_0+0.5,y+0.5,z+0.5  # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
            # print("center node")
            # print(x_c,y_c,z_c)
            graph_yz_0[(x_0_c,y_c,z_c)] = []
            # consider all its neighbors
            for dy, dz in [(0,1),(0,-1),(1,0),(-1,0)]:
                y_, z_ = y_c+dy, z_c+dz
                if 0<y_<n_grid and 0<z_<n_grid:
                    graph_yz_0[(x_0_c,y_c,z_c)].append((x_0_c, y_, z_))
                    # print("center node final")
                    # print("X0YZ")
                    # print("XNYZ")
                    # print(x_0_c,y_,z_)


    y_0=0 #Surface xz at origine
    #We do not redefine shared edges with xy and yz
    for x in range(1, n_grid-1):
        for z in range(1, n_grid-1):
            # print("Node")
            # print(x,y,z)
            x_c,y_0_c,z_c = x+0.5,y_0+0.5,z+0.5 # (x_c,y_c,z_c) is the center of the cube, the real vertices we consider
            # print("center node")
            # print(x_c,y_c,z_c)
            graph_xz_0[(x_c,y_0_c,z_c)] = []
            # consider all its neighbors
            for dx, dz in [(0,1),(0,-1),(1,0),(-1,0)]:
                x_, z_ = x_c+dx, z_c+dz
                if 0<x_<n_grid and 0<z_<n_grid:
                    graph_xz_0[(x_c,y_0_c,z_c)].append((x_, y_0_c, z_))
                    # print("center node final")
                    # print("XY0Z")
                    # print("XYNZ")
                    # print(x_,y_0_c,z_)

    return graph_xy_0, graph_yz_0, graph_xz_0

def get_complete_graph_cube_with_Boundaries(N):
    "For periodic cube-base structure"
    
    G_interior = get_interior_graph_cube(n_grid=N)
    G_exterior_xy_0,G_exterior_yz_0,G_exterior_xz_0  = get_exterior_graph_cube(n_grid=N)
    return G_interior, G_exterior_xy_0,G_exterior_yz_0,G_exterior_xz_0

def get_random_subgraph(graph, p = 0.5):
    """
    Get a random subgraph of a complete graph

    To get the random subgraph, we include each vertex in the complete graph with probability p
    and include all edges in the original graph with two ends selected in the new graph

    input:  graph, the complete graph represented as adjacency list
            p, probability of including each vertex in the subgraph
    output: subgraph, the random subgraph, represented as adjacency list
    """

    vertices_orig =list(graph.keys()) # get all vertices
    vertices = set(random.sample(vertices_orig, int(len(vertices_orig)*p))) # randomly choose with probability p
    # getting the subgraph with the selected vertices
    subgraph = {v:[] for v in vertices}
    for v in vertices:
        for u in graph[v]:
            if u in vertices: # check in u and v both in the new vertex set
                subgraph[v].append(u)
    return subgraph

def get_random_subgraph_pyramid(graph, n_grid, p=0.5):
    """
    Get a random subgraph of a complete graph
    but also add centers of cube surfaces into the graph

    Adding centers make simulation of mechanical properties easier

    input:  graph, the complete graph represented as adjacency list
            p, probability of including each vertex in the subgraph
    output: subgraph, the random subgraph, represented as adjacency list
    """

    # vertices_orig = list(graph.keys())  # get all vertices
    connected_vertices = set()  # set to store vertices that are connected
    main_vertices = set()  # set to store main vertices
    auxiliary_vertices = set() # set to store auxiliary vertices

    # Get all connected vertices in the original graph
    for v, neighbors in graph.items():
        if neighbors:  # if the vertex has neighbors, it's connected
            connected_vertices.add(v)

    for v in connected_vertices:
        coord_x,coord_y,coord_z = v
        if 0 < coord_x < n_grid-1.5 and 0 < coord_y < n_grid-1.5 and 0 < coord_z < n_grid-1.5:
            main_vertices.add(v)
            # print('main')
            # print(v)
        # else:
            # print('auxiliary')
            # print(v)

    # Convert the sets into lists for sampling
    main_vertices_list = list(main_vertices)

    # Randomly select vertices from the main vertices list
    num_selected_main_vertices = int(len(main_vertices_list) * p)
    selected_vertices = random.sample(main_vertices_list, num_selected_main_vertices)

    for v in selected_vertices:
        coord_x,coord_y,coord_z = v
        if coord_x == 0.5:
            auxiliary_vertices.add((coord_x-1,coord_y,coord_z))
        if coord_x == n_grid-2.5:
            auxiliary_vertices.add((coord_x+1,coord_y,coord_z))
        if coord_y == 0.5:
            auxiliary_vertices.add((coord_x,coord_y-1,coord_z))
        if coord_y == n_grid-2.5:
            auxiliary_vertices.add((coord_x,coord_y+1,coord_z))
        if coord_z == 0.5:
            auxiliary_vertices.add((coord_x,coord_y,coord_z-1))
        if coord_z == n_grid-2.5:
            auxiliary_vertices.add((coord_x,coord_y,coord_z+1))

    # Convert the sets into lists for sampling
    auxiliary_vertices_list = list(auxiliary_vertices)

    # Randomly select vertices from the auxiliary vertices list
    num_selected_auxiliary_vertices = int(len(auxiliary_vertices_list) * p)
    selected_vertices.extend(random.sample(auxiliary_vertices_list, num_selected_auxiliary_vertices))

    # Construct the subgraph with the selected vertices
    subgraph = {v: [] for v in selected_vertices}
    for v in selected_vertices:
        for u in graph[v]:
            if u in selected_vertices:
                subgraph[v].append(u)

    return subgraph

def get_random_subgraph_pyramid_1(graph, p=0.5):
    """
    Get a random subgraph of a complete graph
    but also add centers of cube surfaces into the graph

    Adding centers make simulation of mechanical properties easier

    input:  graph, the complete graph represented as adjacency list
            p, probability of including each vertex in the subgraph
    output: subgraph, the random subgraph, represented as adjacency list
    """

    vertices_orig = list(graph.keys())  # get all vertices
    connected_vertices = set()  # set to store vertices that are connected

    # Get all connected vertices in the original graph
    for v, neighbors in graph.items():
        if neighbors:  # if the vertex has neighbors, it's connected
            connected_vertices.add(v)

    # Convert the set of connected vertices into a list for sampling
    connected_vertices_list = list(connected_vertices)

    # Randomly select vertices from the connected vertices list
    num_selected_vertices = int(len(connected_vertices_list) * p)
    selected_vertices = random.sample(connected_vertices_list, num_selected_vertices)

    # Construct the subgraph with the selected vertices
    subgraph = {v: [] for v in selected_vertices}
    for v in selected_vertices:
        for u in graph[v]:
            if u in selected_vertices:
                subgraph[v].append(u)

    return subgraph

def get_random_subgraph_pyramid_old(graph, p = 0.5):
    """
    Get a random subgraph of a complete graph
    but also add centers of cube surfaces into the graph

    Adding centers make simulation of mechanical properties easier

    input:  graph, the complete graph represented as adjacency list
            p, probability of including each vertex in the subgraph
    output: subgraph, the random subgraph, represented as adjacency list
    """

    vertices_orig =list(graph.keys())  # get all vertices
    #print("Existing keys:", graph.keys())
    vertices = set(random.sample(vertices_orig, int(len(vertices_orig)*p))) # randomly choose with probability p

    # getting the subgraph with the selected vertices
    subgraph = {v:[] for v in vertices}
    for v in vertices:
        for u in graph[v]:
            if u in vertices:
                subgraph[v].append(u)
    return subgraph

def get_random_spanning_tree(graph):
    """
    Get a random spanning tree from a graph

    A spanning tree is a spanning tree T of an undirected graph G is a subgraph that is a tree which includes all of the vertices of G
    This is just one way of implementing this, see the following link for more discussions
    https://stackoverflow.com/questions/2041517/random-simple-connected-graph-generation-with-given-sparseness


    input:  graph, the orginal graph, represented as adjacency list
    output: span_tree, the spanning tree, represented as adjacency list
            return None if the original graph is not connected
    """


    vertices = set(graph.keys()) # all vertices in the graph
    v0 = vertices.pop() # get the first vertex
    span_tree = {v0:[]} # initialize the spanning tree
    # the neighbors_dict and neighbors are used to maintain a random set structure
    # see https://leetcode.com/problems/insert-delete-getrandom-o1/
    neighbors_dict = {v:i for i,v in enumerate(list(graph[v0]))}
    neighbors = list(graph[v0])

    while neighbors:
        v = random.choice(neighbors) # randomly choose a neighbor
        # the following 5 lines try to delete v from the neighbor set with O(1) complexity
        i = neighbors_dict[v]
        neighbors[i], neighbors[-1] = neighbors[-1], neighbors[i]
        neighbors_dict[neighbors[i]] = i
        neighbors.pop()
        neighbors_dict.pop(v)

        vertices.remove(v) # also remove v from vertices
        span_tree[v] = []
        v_neighbors = list(graph[v]) # get all neighbors of v
        random.shuffle(v_neighbors) # random shuffle
        for u in v_neighbors:
            if u in span_tree:
                # add the edge to spanning tree, the first u in neighbors of v that is include in the spanning tree
                if not span_tree[v]:
                    span_tree[v].append(u)
                    span_tree[u].append(v)
            else:
                # include vertices adjacent to v into the neighbor set
                if u not in neighbors_dict:
                    neighbors.append(u)
                    neighbors_dict[u] = len(neighbors)-1

    # in the vertices is not empty, it means the input graph is not connected
    if vertices:
        # print(vertices)
        print("The input graph is not connected")
        return None
    else:
        return span_tree

def add_random_complement_graph(graph_sub, graph_orig, p = 0.5):
    """
    Adding random edges to a subgraph

    edges in the original graph but not in the subgraph are added with probability p

    input:  graph_sub, the subgraph, represented as adjacency list
            graph_orig, the original graph, represented as adjacency list
    output: graph_added, the new graph, represented as adjacency list
    """

    graph_added = graph_sub.copy()
    edges_complement = [] # the set of edges not in subgraph but in the original graph
    for v in graph_orig:
        for u in graph_orig[v]:
            if v<u and u not in graph_sub[v]:
                # only consider (v,u) but not (u,v) since its undirected graph
                # also here we check if (v,u) is in the subgraph, if not, we add this into the graph
                edges_complement.append((v,u))
    edges_add = random.sample(edges_complement, int(p*len(edges_complement))) # randomly select the complement edges with probability p
    # for edge in edges_add:
    #     print("End nodes of edge:", edge)

    # adding the select edges
    for (v,u) in edges_add:
        graph_added[v].append(u)
        graph_added[u].append(v)

    return graph_added


def add_aligned_edges(G_span_tree, G_randsub, p, bx, by, bz):
    """
    Adding edges aligned in the desiered direction to a subgraph (spanning tree)
    Edges in the original graph but not in the subgraph and aligned in the z-direction are added.
    input:  G_span_tree, the subgraph (spanning tree), represented as adjacency list
            G_randsub, the original graph (first subgraph), represented as adjacency list
    output: graph_added, the new graph, represented as adjacency list
    """

    graph_added = G_span_tree.copy()

    edges_aligned_x = []  # the set of edges not in subgraph but in the original graph and aligned in x-direction
    edges_aligned_y = []
    edges_aligned_z = []

    edges_aligned_xy = []
    edges_aligned_yz = []
    edges_aligned_xz = []

    # Iterate over the original graph to find complement edges
    for v in G_randsub.keys():
        for u in G_randsub[v]:
            # Check if the edge (v, u) is not already in the subgraph and v < u to avoid duplicate edges
            if v < u and u not in G_span_tree.get(v, []):
                # Check if the difference in x-coordinates is 1 and differences in other dimensions are 0

                if abs(v[0] - u[0]) == 1 and abs(v[1] - u[1]) == 0 and abs(v[2] - u[2]) == 0:
                    edges_aligned_x.append((v, u))
                if abs(v[1] - u[1]) == 1 and abs(v[0] - u[0]) == 0 and abs(v[2] - u[2]) == 0:
                    edges_aligned_y.append((v, u))
                if abs(v[2] - u[2]) == 1 and abs(v[0] - u[0]) == 0 and abs(v[1] - u[1]) == 0:
                    edges_aligned_z.append((v, u))

                if (abs(v[0] - u[0]) == 1 or abs(v[1] - u[1]) == 1) and abs(v[2] - u[2]) == 0:
                    edges_aligned_xy.append((v, u))
                if (abs(v[1] - u[1]) == 1 or abs(v[2] - u[2]) == 1) and abs(v[0] - u[0]) == 0:
                    edges_aligned_yz.append((v, u))
                if (abs(v[0] - u[0]) == 1 or abs(v[2] - u[2]) == 1) and abs(v[1] - u[1]) == 0:
                    edges_aligned_xz.append((v, u))

    edges_add_x = random.sample(edges_aligned_x, int(p*bx*len(edges_aligned_x)))
    edges_add_y = random.sample(edges_aligned_y, int(p*by*len(edges_aligned_y)))
    edges_add_z = random.sample(edges_aligned_z, int(p*bz*len(edges_aligned_z)))

    edges_add_xy = random.sample(edges_aligned_xy, int(p*bx*len(edges_aligned_xy)))
    edges_add_yz = random.sample(edges_aligned_yz, int(p*by*len(edges_aligned_yz)))
    edges_add_xz = random.sample(edges_aligned_xz, int(p*bz*len(edges_aligned_xz)))

    for (v,u) in edges_add_x:
        graph_added[v].append(u)
        graph_added[u].append(v)

    for (v,u) in edges_add_y:
        graph_added[v].append(u)
        graph_added[u].append(v)

    for (v,u) in edges_add_z:
        graph_added[v].append(u)
        graph_added[u].append(v)

    for (v,u) in edges_add_xy:
        graph_added[v].append(u)
        graph_added[u].append(v)

    for (v,u) in edges_add_yz:
        graph_added[v].append(u)
        graph_added[u].append(v)

    for (v,u) in edges_add_xz:
        graph_added[v].append(u)
        graph_added[u].append(v)

    return graph_added

def get_rods_from_graph(graph):
    """
    Calculate the edges of the graph
    input: the structure represented as a graph
    output: ls, a list of all edges in the graph
    """
    ls = []
    for p1 in graph.keys():
        for p2 in graph[p1]:
            if p1 < p2:
                ls.append(list(p1) + list(p2))
    ls = np.array(ls)
    ls = ls.reshape((-1, 2, 3))
    # print('number of edges:', len(ls))
    return ls

def get_triangles_from_graph(graph):
    """
    Get the triangles from plane-based structures
    input:  graph, represented as adjacency list
    output: triangles, list of 9 element array containing the xyz coordinates of the vertices of the triangle
    """
    triangles = set() # it's a set because we might add reduandant triangles, this will remove redundancy
    for p1 in graph.keys():
        for p2 in graph[p1]:
            if p1<p2:
                x1, y1, z1 = p1
                x2, y2, z2 = p2
                n_diff = int(x1!=x2) + int(y1!=y2) + int(z1!=z2)
                if n_diff==1:
                    if x1!=x2:
                        if int(y1)==y1:
                            x3, y3, z3 = (x1+x2)/2, y1, z1-0.5
                            x4, y4, z4 = (x1+x2)/2, y1, z1+0.5
                        else:
                            x3, y3, z3 = (x1+x2)/2, y1-0.5, z1
                            x4, y4, z4 = (x1+x2)/2, y1+0.5, z1
                    elif y1!=y2:
                        if int(x1)==x1:
                            x3, y3, z3 = x1, (y1+y2)/2, z1-0.5
                            x4, y4, z4 = x1, (y1+y2)/2, z1+0.5
                        else:
                            x3, y3, z3 = x1-0.5, (y1+y2)/2, z1
                            x4, y4, z4 = x1+0.5, (y1+y2)/2, z1
                    elif z1!=z2:
                        if int(x1)==x1:
                            x3, y3, z3 = x1, y1-0.5, (z1+z2)/2
                            x4, y4, z4 = x1, y1+0.5, (z1+z2)/2
                        else:
                            x3, y3, z3 = x1-0.5, y1, (z1+z2)/2
                            x4, y4, z4 = x1+0.5, y1, (z1+z2)/2
                elif n_diff==2:
                    if x1==x2:
                        if int(y1)==y1:
                            x3, y3, z3 = x1-0.5, y1, z2
                            x4, y4, z4 = x1+0.5, y1, z2
                        else:
                            x3, y3, z3 = x1-0.5, y2, z1
                            x4, y4, z4 = x1+0.5, y2, z1
                    elif y1==y2:
                        if int(x1)==x1:
                            x3, y3, z3 = x1, y1-0.5, z2
                            x4, y4, z4 = x1, y1+0.5, z2
                        else:
                            x3, y3, z3 = x2, y1-0.5, z1
                            x4, y4, z4 = x2, y1+0.5, z1
                    elif z1==z2:
                        if int(x1)==x1:
                            x3, y3, z3 = x1, y2, z1-0.5
                            x4, y4, z4 = x1, y2, z1+0.5
                        else:
                            x3, y3, z3 = x2, y1, z1-0.5
                            x4, y4, z4 = x2, y1, z1+0.5
                else:
                    assert(0)
                # (x3, y3, z3) and (x4, y4, z4) are the coordicates of the shared edge of two squares
                triangles.add((x1, y1, z1, x3, y3, z3, x4, y4, z4))
                triangles.add((x2, y2, z2, x3, y3, z3, x4, y4, z4))

    triangles = np.array(list(triangles))

    return triangles

def get_pyramids_from_graph(graph, Added_Grid_Number):
    """
    Get the pyramids from pyramid-based structures
    input:  graph, represented as adjacency list
    output: pyramids, list of 15 element array containing the xyz coordinates of the vertices of the pyramid
    """
    pyramids = set() # it's a set because we might add reduandant triangles, this will remove redundancy
    for p1 in graph.keys():
        for p2 in graph[p1]:
            if p1<p2:
                x1, y1, z1 = p1
                x2, y2, z2 = p2
                #print(p1)
                #print(p2)
                #calculating the corners for the pyramids and removing the pyramids over the range
                if all(0 < vv < Added_Grid_Number-1.5 for vv in [x1, y1, z1]):
                    #print("One")
                    if x1!=x2:
                        #print("Two")
                        x3, y3, z3 = (x1+x2)/2, y1+0.5, z1+0.5
                        x4, y4, z4 = (x1+x2)/2, y1-0.5, z1+0.5
                        x5, y5, z5 = (x1+x2)/2, y1-0.5, z1-0.5
                        x6, y6, z6 = (x1+x2)/2, y1+0.5, z1-0.5
                    elif y1!=y2:
                        #print("Three")
                        x3, y3, z3 = x1-0.5, (y1+y2)/2, z1+0.5
                        x4, y4, z4 = x1+0.5, (y1+y2)/2, z1+0.5
                        x5, y5, z5 = x1+0.5, (y1+y2)/2, z1-0.5
                        x6, y6, z6 = x1-0.5, (y1+y2)/2, z1-0.5
                    elif z1!=z2:
                        #print("Four")
                        x3, y3, z3 = x1-0.5, y1-0.5, (z1+z2)/2
                        x4, y4, z4 = x1+0.5, y1-0.5, (z1+z2)/2
                        x5, y5, z5 = x1+0.5, y1+0.5, (z1+z2)/2
                        x6, y6, z6 = x1-0.5, y1+0.5, (z1+z2)/2
                    else:
                        #print("Five")
                        assert(0)
                    # (x3, y3, z3), (x4, y4, z4), (x5, y, z5), (x6, y6, z6) are the coordicates of the shared surface of two pyramids
                    pyramids.add((x1, y1, z1, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6))
                    #pyramids.add((x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6))
                if all(0 < vv < Added_Grid_Number-1.5 for vv in [x2, y2, z2]):
                    #print("One")
                    if x1!=x2:
                        #print("Two")
                        x3, y3, z3 = (x1+x2)/2, y1+0.5, z1+0.5
                        x4, y4, z4 = (x1+x2)/2, y1-0.5, z1+0.5
                        x5, y5, z5 = (x1+x2)/2, y1-0.5, z1-0.5
                        x6, y6, z6 = (x1+x2)/2, y1+0.5, z1-0.5
                    elif y1!=y2:
                        #print("Three")
                        x3, y3, z3 = x1-0.5, (y1+y2)/2, z1+0.5
                        x4, y4, z4 = x1+0.5, (y1+y2)/2, z1+0.5
                        x5, y5, z5 = x1+0.5, (y1+y2)/2, z1-0.5
                        x6, y6, z6 = x1-0.5, (y1+y2)/2, z1-0.5
                    elif z1!=z2:
                        #print("Four")
                        x3, y3, z3 = x1-0.5, y1-0.5, (z1+z2)/2
                        x4, y4, z4 = x1+0.5, y1-0.5, (z1+z2)/2
                        x5, y5, z5 = x1+0.5, y1+0.5, (z1+z2)/2
                        x6, y6, z6 = x1-0.5, y1+0.5, (z1+z2)/2
                    else:
                        #print("Five")
                        assert(0)
                    # (x3, y3, z3), (x4, y4, z4), (x5, y, z5), (x6, y6, z6) are the coordicates of the shared surface of two pyramids
                    pyramids.add((x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6))

    # # adding two the entire surfaces of two boundaries in Z direction
    # for x in range(Added_Grid_Number-2):
    #     for y in range(Added_Grid_Number-2):
    #         #print((x+0.5, y+0.5, 0.5))
    #         #print((x+0.5, y+0.5, N-1.5))
    #         x1, y1 = x+0.5, y+0.5
    #         x3, y3 = x1-0.5, y1-0.5
    #         x4, y4 = x1+0.5, y1-0.5
    #         x5, y5 = x1+0.5, y1+0.5
    #         x6, y6 = x1-0.5, y1+0.5
    #         #print((x1, y1, 0.5, x3, y3, 0, x4, y4, 0, x5, y5, 0, x6, y6, 0))
    #         #print((x1, y1, N-2.5, x3, y3, N-2, x4, y4, N-2, x5, y5, N-2, x6, y6, N-2))
    #         pyramids.add((x1, y1, 0.5, x3, y3, 0, x4, y4, 0, x5, y5, 0, x6, y6, 0))
    #         pyramids.add((x1, y1, Added_Grid_Number-2.5, x3, y3, Added_Grid_Number-2, x4, y4, Added_Grid_Number-2, x5, y5, Added_Grid_Number-2, x6, y6, Added_Grid_Number-2))

    pyramids = np.array(list(pyramids))

    return pyramids

def get_cubes_from_graph(graph, Grid_Number):
    """
    Get the cubes from cube-based structures
    input:  graph, represented as adjacency list
    output: cubes, list of 3 element array containing the xyz coordinates of the center node of a cube
    """
    cubes = set() # it's a set because we might add redundant nodes, this will remove redundancy
    for p1 in graph.keys():
        for p2 in graph[p1]:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            cubes.add((x1, y1, z1))
            cubes.add((x2, y2, z2))
    # #Adding the whole face in Z direction
    # for x in range(Grid_Number):
    #     for y in range(Grid_Number):
    #         cubes.add((x+0.5, y+0.5, 0.5))
    #         cubes.add((x+0.5, y+0.5, Grid_Number-0.5))

    cubes = np.array(list(cubes))
    return cubes

def get_cubes_from_graph_with_boundaries(N, G_interior_randadded, G_exterior_xy_0_randadded, G_exterior_yz_0_randadded, G_exterior_xz_0_randadded):
    """
    Get the cubes from cube-based structures
    input:  graph, represented as adjacency list
    output: cubes, list of 3 element array containing the xyz coordinates of the center node of a cube
    """
    cubes = set() # it's a set because we might add redundant nodes, this will remove redundancy
    for p1 in G_interior_randadded.keys():
        for p2 in G_interior_randadded[p1]:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            cubes.add((x1, y1, z1))
            cubes.add((x2, y2, z2))

    for p1 in G_exterior_xy_0_randadded.keys():
        for p2 in G_exterior_xy_0_randadded[p1]:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            cubes.add((x1, y1, z1))
            cubes.add((x2, y2, z2))
            cubes.add((x1, y1, z1+N-1)) #for the opposite surface
            cubes.add((x2, y2, z2+N-1)) #for the opposite surface

    for p1 in G_exterior_yz_0_randadded.keys():
        for p2 in G_exterior_yz_0_randadded[p1]:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            cubes.add((x1, y1, z1))
            cubes.add((x2, y2, z2))
            cubes.add((x1+N-1, y1, z1)) #for the opposite surface
            cubes.add((x2+N-1, y2, z2)) #for the opposite surface

    for p1 in G_exterior_xz_0_randadded.keys():
        for p2 in G_exterior_xz_0_randadded[p1]:
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            cubes.add((x1, y1, z1))
            cubes.add((x2, y2, z2))
            cubes.add((x1, y1+N-1, z1)) #for the opposite surface
            cubes.add((x2, y2+N-1, z2)) #for the opposite surface

    #adding edges in full
    for x in range(N):
        cubes.add((x+0.5, 0.5, 0.5))
        cubes.add((x+0.5, N-0.5, 0.5))
        cubes.add((x+0.5, 0.5, N-0.5))
        cubes.add((x+0.5, N-0.5, N-0.5))

    for y in range(N):
        cubes.add((0.5, y+0.5, 0.5))
        cubes.add((N-0.5, y+0.5, 0.5))
        cubes.add((0.5, y+0.5, N-0.5))
        cubes.add((N-0.5, y+0.5, N-0.5))

    for z in range(N):
        cubes.add((0.5, 0.5, z+0.5))
        cubes.add((N-0.5, N-0.5, z+0.5))
        cubes.add((0.5, N-0.5, z+0.5))
        cubes.add((N-0.5, 0.5, z+0.5))

    cubes = np.array(list(cubes))
    return cubes

def plot_3d_graph(graph):
    """
    Plot the rod-based structures
    input: the structure represented as a graph
    output: ls, a list of all edges in the graph
            + we also plot the edges
    """
    fig = plt.figure()
    # ax = fig.gca(projection='3d')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$')
    ax.view_init(elev=20., azim=20)
    pls = np.array(list(graph.keys()))
    ax.plot(pls[:,0], pls[:,1], pls[:,2], '.r', markersize=3)
    ls = []
    for p1 in graph.keys():
        for p2 in graph[p1]:
            if p1<p2:
                ls.append(list(p1)+list(p2))
    ls = np.array(ls)
    # print('number of edges:', len(ls))
    ls = ls.reshape((-1,2,3))
    lc = Line3DCollection(ls, linewidths=0.5, colors='b')
    ax.add_collection(lc)
    plt.show()
    return ls

def plot_3d_graph_noName(graph, geometry, grid_number, filename):
    """
    Plot the rod-based structures
    input: the structure represented as a graph
    output: ls, a list of all edges in the graph
            + we also plot the edges, a cube, and grids
    """
    fig = plt.figure(figsize=(5, 8)) # Set figure size to 3x3 inches
    ax = fig.add_subplot(111, projection='3d')

    # Remove the borders, axis labels, ticks, and numbers
    ax.set_axis_off()

    # Set the view angle
    ax.view_init(elev=25, azim=230)
    ax.set_box_aspect([1, 1, 1.01])  # Aspect ratio is x:y:z


    # Calculate the cube corner at (0.5, 0.5, 0.5)
    # cube_corner = np.array([0,0,0])

    # Generate grid points
    if geometry == 'Pyramid':
        grid_points = np.arange(-1, grid_number + 2)
    else:
        grid_points = np.arange(0, grid_number)
    # grid_points = np.arange(0, grid_number)
    # print(grid_points)

    # Plot the grids with light gray dashed lines
    if geometry == 'Pyramid':
        for x in grid_points:
            for y in grid_points:
                ax.plot3D([x, x], [y, y], [-1, grid_number+1], color="lightgray", linestyle="--", linewidth=0.5)
                ax.plot3D([x, x], [-1, grid_number+1], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
                ax.plot3D([-1, grid_number+1], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
    else:
        for x in grid_points:
            for y in grid_points:
                ax.plot3D([x, x], [y, y], [0, grid_number-1], color="lightgray", linestyle="--", linewidth=0.5)
                ax.plot3D([x, x], [0, grid_number-1], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
                ax.plot3D([0, grid_number-1], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)

    # Generate grid points
    # grid_points = np.linspace(0, (grid_number-1), grid_number)
    # grid_points = np.linspace(0, cube_size, grid_number)
    #
    # # Plot the grids with light gray dashed lines
    # for x in grid_points:
    #     for y in grid_points:
    #         ax.plot3D([x, x], [y, y], [0, cube_size], color="lightgray", linestyle="--", linewidth=0.5)
    #         ax.plot3D([x, x], [0, cube_size], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
    #         ax.plot3D([0, cube_size], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
    #         # ax.plot3D([x, x], [y, y], [0, (grid_number-1)], color="lightgray", linestyle="--", linewidth=0.5)
    #         # ax.plot3D([x, x], [0, (grid_number-1)], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
    #         # ax.plot3D([0, (grid_number-1)], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
    #
    # # Calculate the cube size based on the maximum range of coordinates
    # max_range = max(max(p) for p in graph.keys())
    # cube_size = max_range / grid_number
    #
    # # # Plot the cube
    # # for s, e in [(s, e) for s in [0, cube_size] for e in [0, cube_size] if s != e]:
    # #     ax.plot3D([s, e], [0, 0], [0, 0], color="lightgray", linestyle="--", linewidth=0.5)
    # #     ax.plot3D([0, 0], [s, e], [0, 0], color="lightgray", linestyle="--", linewidth=0.5)
    # #     ax.plot3D([0, 0], [0, 0], [s, e], color="lightgray", linestyle="--", linewidth=0.5)
    # #
    # #     ax.plot3D([s, e], [0, 0], [cube_size, cube_size], color="lightgray", linestyle="--", linewidth=0.5)
    # #     ax.plot3D([0, 0], [s, e], [cube_size, cube_size], color="lightgray", linestyle="--", linewidth=0.5)
    # #     ax.plot3D([0, 0], [0, 0], [s, e], color="lightgray", linestyle="--", linewidth=0.5)
    # #
    # #     ax.plot3D([0, 0], [s, e], [0, 0], color="lightgray", linestyle="--", linewidth=0.5)
    # #     ax.plot3D([s, e], [0, 0], [0, 0], color="lightgray", linestyle="--", linewidth=0.5)
    #
    # print("Graph Points Coordinates:")
    # for p in graph.keys():
    #     print(f"Point: {p}")

    # Plot the nodes and edges
    for p1 in graph.keys():
        for p2 in graph[p1]:
            ax.plot3D([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='b', linewidth=0.5)
    for p in graph.keys():
        ax.scatter(p[0], p[1], p[2], color='r', s=1)  # Adjust marker size here

    # plt.show()
    plt.savefig(filename, dpi=1000, bbox_inches='tight', format='png')
    plt.close()


# def plot_3d_graph_noNameCube(graph, cube_size, grid_number, filename):
#     """
#     Plot the rod-based structures
#     input: the structure represented as a graph
#     output: ls, a list of all edges in the graph
#             + we also plot the edges, a cube, and grids
#     """
#     fig = plt.figure(figsize=(5, 8))
#     ax = fig.add_subplot(111, projection='3d')
#
#     # Remove the borders, axis labels, ticks, and numbers
#     ax.set_axis_off()
#
#     # Set the view angle
#     ax.view_init(elev=25, azim=230)
#     ax.set_box_aspect([1, 1, 1.01])  # Aspect ratio is x:y:z
#
#     # Calculate the cube corner at (0.5, 0.5, 0.5)
#     cube_corner = np.array([0,0,0])
#
#     # Generate grid points
#     grid_points = np.arange(0, grid_number)
#     # print(grid_points)
#
#     # # Plot the grids with light gray dashed lines
#     # for x in grid_points:
#     #     for y in grid_points:
#     #         ax.plot3D([x, x], [y, y], [cube_corner[2], cube_corner[2] + cube_size], color="lightgray", linestyle="--", linewidth=0.5)
#     #         ax.plot3D([x, x], [cube_corner[1], cube_corner[1] + cube_size], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
#     #         ax.plot3D([cube_corner[0], cube_corner[0] + cube_size], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
#     #
#     # # Plot the cube
#     # for s, e in [(s, e) for s in [cube_corner[0], cube_corner[0] + cube_size] for e in [cube_corner[1], cube_corner[1] + cube_size] if s != e]:
#     #     ax.plot3D([s, e], [cube_corner[1], cube_corner[1]], [cube_corner[2], cube_corner[2]], color="lightgray", linestyle="--", linewidth=0.5)
#     #     ax.plot3D([cube_corner[0], cube_corner[0]], [s, e], [cube_corner[2], cube_corner[2]], color="lightgray", linestyle="--", linewidth=0.5)
#     #     ax.plot3D([cube_corner[0], cube_corner[0]], [cube_corner[1], cube_corner[1]], [s, e], color="lightgray", linestyle="--", linewidth=0.5)
#     #
#     #     ax.plot3D([s, e], [cube_corner[1], cube_corner[1]], [cube_corner[2] + cube_size, cube_corner[2] + cube_size], color="lightgray", linestyle="--", linewidth=0.5)
#     #     ax.plot3D([cube_corner[0], cube_corner[0]], [s, e], [cube_corner[2] + cube_size, cube_corner[2] + cube_size], color="lightgray", linestyle="--", linewidth=0.5)
#     #     ax.plot3D([cube_corner[0], cube_corner[0]], [cube_corner[1], cube_corner[1]], [s, e], color="lightgray", linestyle="--", linewidth=0.5)
#     #
#     #     ax.plot3D([cube_corner[0], cube_corner[0]], [s, e], [cube_corner[2], cube_corner[2]], color="lightgray", linestyle="--", linewidth=0.5)
#     #     ax.plot3D([s, e], [cube_corner[1], cube_corner[1]], [cube_corner[2], cube_corner[2]], color="lightgray", linestyle="--", linewidth=0.5)
#         # Plot the grids with light gray dashed lines
#         # Generate grid points
#
#     # Plot the grids with light gray dashed lines
#     for x in grid_points:
#         for y in grid_points:
#             ax.plot3D([x, x], [y, y], [0, grid_number-1], color="lightgray", linestyle="--", linewidth=0.5)
#             ax.plot3D([x, x], [0, grid_number-1], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
#             ax.plot3D([0, grid_number-1], [x, x], [y, y], color="lightgray", linestyle="--", linewidth=0.5)
#
#     # Plot the nodes and edges
#     for p1 in graph.keys():
#         for p2 in graph[p1]:
#             ax.plot3D([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='b', linewidth=0.5)
#     for p in graph.keys():
#         ax.scatter(p[0], p[1], p[2], color='r', s=1)  # Adjust marker size here
#
#     print("Graph Points Coordinates:")
#     for p in graph.keys():
#         print(f"Point: {p}")
#
#     # plt.show()
#     plt.savefig(filename, dpi=1000, bbox_inches='tight', format='png')
#     plt.close()

def plot_3d_2D_graph(graph):
    """
    Plot the rod-based structures
    input: the structure represented as a graph
    output: ls, a list of all edges in the graph
            + we also plot the edges
    """
    fig = plt.figure(figsize=(10, 5))

    # 3D Plot
    ax_3d = fig.add_subplot(121, projection='3d')
    ax_3d.set_xlabel('$X$')
    ax_3d.set_ylabel('$Y$')
    ax_3d.set_zlabel('$Z$')

    pls = np.array(list(graph.keys()))
    ax_3d.plot(pls[:,0], pls[:,1], pls[:,2], '.r', markersize=3)

    ls = []
    for p1 in graph.keys():
        for p2 in graph[p1]:
            if p1 < p2:
                ls.append(list(p1)+list(p2))
    ls = np.array(ls)
    # print('number of edges:', len(ls))
    ls = ls.reshape((-1, 2, 3))
    lc = Line3DCollection(ls, linewidths=0.5, colors='b')
    ax_3d.add_collection(lc)

    # XY Plane Plot
    ax_xy = fig.add_subplot(122)
    ax_xy.set_xlabel('$X$')
    ax_xy.set_ylabel('$Y$')

    ax_xy.plot(pls[:,0], pls[:,1], '.r', markersize=3)
    for edge in ls:
        ax_xy.plot(edge[:,0], edge[:,1], color='b', linewidth=0.5)

    plt.tight_layout()
    plt.show()

    return ls
##-----------------------Run Functions-----------------------------------------##
def Rod_Run(p1, p2, bx, by, bz, Grid_Number, cube_size, NameData):
  "generate rod based random structure"

  # Set the output file name
  filename = NameData+".mat"
  # filename_z = NameData+"_z.mat"
  # full_path = path + "\\" + filename
  full_path = os.path.join(path, filename)
  # full_path_z = os.path.join(path, filename_z)

  G_complete = get_complete_graph_rod(n_grid=Grid_Number)
  # plot_3d_graph(G_complete)
  plot_3d_graph_noName(G_complete, 'Rod', Grid_Number, 'RodG_complete.png')

  # plot_3d_2D_graph(G_complete)
  # print(lss)
  # ls_C = get_rods_from_graph(G_complete)
  G_randsub = get_random_subgraph(G_complete,p=p1)
  # plot_3d_graph(G_randsub)
  plot_3d_graph_noName(G_randsub, 'Rod', Grid_Number, 'RodG_randsub.png')

  # ls_G = get_rods_from_graph(G_randsub)
  G_span_tree = get_random_spanning_tree(G_randsub)
  # ls_S = get_rods_from_graph(G_span_tree)
  # plot_3d_graph(G_span_tree)
  plot_3d_graph_noName(G_span_tree, 'Rod', Grid_Number, 'RodG_span_tree.png')


  # G_randadded = add_random_complement_graph(G_span_tree, G_randsub, p2)
  G_randadded = add_aligned_edges(G_span_tree, G_randsub, p2, bx, by, bz)
  # plot_3d_2D_graph(G_randadded)
  # plot_3d_graph(G_randadded)
  plot_3d_graph_noName(G_randadded, 'Rod', Grid_Number, 'RodG_randadded.png')

  ls_b = get_rods_from_graph(G_randadded)
  sio.savemat(full_path, {'edges': ls_b*cube_size/(Grid_Number-1), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})
  # sio.savemat(full_path, {'edges': ls_b*cube_size/(Grid_Number-1), 'p1': p1, 'p2': p2, 'Grid_Number': Grid_Number})

  # G_randaddedx = add_aligned_x_edges(G_span_tree,G_randsub,p=p2)
  # ls_x = get_rods_from_graph(G_randaddedx)
  # sio.savemat(full_path, {'edges': ls_x*cube_size/(Grid_Number-1)})

  # G_randaddedy = add_aligned_y_edges(G_span_tree,G_randsub,p=p2)
  # ls_y = get_rods_from_graph(G_randaddedy)
  # sio.savemat(full_path, {'edges': ls_y*cube_size/(Grid_Number-1)})

  # G_randaddedz = add_aligned_z_edges(G_span_tree,G_randsub,p=p2)
  # ls_z = get_rods_from_graph(G_randaddedz)
  # sio.savemat(full_path, {'edges': ls_z*cube_size/(Grid_Number-1)})

  # G_randaddedxy = add_aligned_xy_edges(G_span_tree,G_randsub,p=p2)
  # ls_xy = get_rods_from_graph(G_randaddedxy)
  # sio.savemat(full_path, {'edges': ls_xy*cube_size/(Grid_Number-1)})

  # G_randaddedxz = add_aligned_xz_edges(G_span_tree,G_randsub,p=p2)
  # ls_xz = get_rods_from_graph(G_randaddedxz)
  # sio.savemat(full_path, {'edges': ls_xz*cube_size/(Grid_Number-1)})

  # G_randaddedyz = add_aligned_yz_edges(G_span_tree,G_randsub,p=p2)
  # ls_yz = get_rods_from_graph(G_randaddedyz)
  # sio.savemat(full_path, {'edges': ls_yz*cube_size/(Grid_Number-1)})

  # G_randadded = add_random_complement_graph(G_span_tree,G_randsub,p=p2)
  # ls = get_rods_from_graph(G_randadded)
  # sio.savemat(full_path, {'edges': ls*cube_size/(Grid_Number-1)})

  pass

def Triangle_Run(p1, p2, bx, by, bz, Grid_Number, cube_size, NameData):
  "generate triangular based random structure"

  # Set the output file name
  filename = NameData+".mat"
  # full_path = path + "\\" + filename
  full_path = os.path.join(path, filename)

  G_complete = get_complete_graph_plane(n_grid=Grid_Number)
  # plot_3d_2D_graph(G_complete)
  plot_3d_graph_noName(G_complete, 'Triangle', Grid_Number,'TriangleG_complete.png')
  G_randsub = get_random_subgraph(G_complete,p=p1)
  plot_3d_graph_noName(G_randsub, 'Triangle', Grid_Number,'TriangleG_randsub.png')
  G_span_tree = get_random_spanning_tree(G_randsub)
  plot_3d_graph_noName(G_span_tree, 'Triangle', Grid_Number,'TriangleG_span_tree.png')
  G_randadded = add_random_complement_graph(G_span_tree,G_randsub,p=p2)
  plot_3d_graph_noName(G_randadded, 'Triangle', Grid_Number,'TriangleG_randadded.png')
  # G_randadded = add_aligned_edges(G_span_tree, G_randsub, p2, bx, by, bz)
  triangles = get_triangles_from_graph(G_randadded)
  sio.savemat(full_path, {'triangles': triangles*cube_size/(Grid_Number-1), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})
  pass

def Pyramid_Run(p1, p2, bx, by, bz, Grid_Number, cube_size, NameData):
  "generate pyramid/truncated pyramid based random structure"
  Grid_Number = Grid_Number -1
  # Set the output file name
  filename = NameData+".mat"
  # full_path = path + "\\" + filename
  full_path = os.path.join(path, filename)

  Added_Grid_Number = Grid_Number + 2 # We include 2 more nodes in each direction to include boundary pyramids as well. We remove them at the end

  G_complete = get_complete_graph_pyramid_old(n_grid=Added_Grid_Number)
  # G_complete = get_complete_graph_pyramid(n_grid=Added_Grid_Number)
  pyramidsC = get_pyramids_from_graph(G_complete, Added_Grid_Number)
  sio.savemat(full_path, {'pyramids': pyramidsC*cube_size/(Grid_Number), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})

  # plot_3d_2D_graph(G_complete)
  plot_3d_graph_noName(G_complete, 'Pyramid', Grid_Number,'PyramidG_complete.png')
  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=20:
      # G_randsub = get_random_subgraph_pyramid(G_complete, Added_Grid_Number, p=p1)
      G_randsub = get_random_subgraph_pyramid_old(G_complete, p=p1)
      G_span_tree = get_random_spanning_tree(G_randsub)
      if G_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1

  plot_3d_graph_noName(G_randsub, 'Pyramid', Grid_Number, 'PyramidG_randsub.png')
  # plot_3d_2D_graph(G_span_tree)
  plot_3d_graph_noName(G_span_tree, 'Pyramid', Grid_Number, 'PyramidG_span_tree.png')
  G_randadded = add_random_complement_graph(G_span_tree,G_randsub,p=p2)
  G_randadded = add_aligned_edges(G_span_tree, G_randsub, p2, bx, by, bz)
  # plot_3d_2D_graph(G_randadded)
  plot_3d_graph_noName(G_randadded, 'Pyramid', Grid_Number, 'PyramidG_randadded.png')
  # print('number of edges', ls.shape)
  pyramids = get_pyramids_from_graph(G_randadded, Added_Grid_Number)
  sio.savemat(full_path, {'pyramids': pyramids*cube_size/(Grid_Number), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})

  #eng = matlab.engine.start_matlab()
  #connect = eng.comsol_start()

  #Impedance, VolumeFraction = eng.Pyramid_Solid(filename, cube_size)
  # Impedance = eng.Pyramid_Solid(filename, cube_size)
  # Results = eng.Pyramid_Solid(filename, cube_size)

  #sio.savemat(full_path, {'impedance': Impedance})

  # print('number of pyramids', pyramids.shape)
  # VolumeFraction = pyramids.shape * 1/3 * 1/2 * (cube_size/(Grid_Number))^3
  # sio.savemat(full_path, {'volumefraction': VolumeFraction})

  pass

def Cube_Run(p1, p2, bx, by, bz, Grid_Number, cube_size, NameData):
  "generate cube based random structure"
  Grid_Number = Grid_Number -1
  # Set the output file name
  filename = NameData+".mat"
  # full_path = path + "\\" + filename
  full_path = os.path.join(path, filename)

  G_complete = get_complete_graph_cube(n_grid=Grid_Number)
  # plot_3d_2D_graph(G_complete)
  plot_3d_graph_noName(G_complete, 'Cube', Grid_Number+1,'CubeG_complete.png')
  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=10:
      G_randsub = get_random_subgraph(G_complete,p=p1)
      G_span_tree = get_random_spanning_tree(G_randsub)
      if G_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1
  G_randadded = add_random_complement_graph(G_span_tree,G_randsub,p=p2)
  plot_3d_graph_noName(G_randsub, 'Cube', Grid_Number+1, 'CubeG_randsub.png')
  plot_3d_graph_noName(G_span_tree, 'Cube', Grid_Number+1, 'CubeG_span_tree.png')
  # G_randadded = add_aligned_edges(G_span_tree, G_randsub, p2, bx, by, bz)
  plot_3d_graph_noName(G_randadded, 'Cube', Grid_Number+1, 'CubeG_randadded.png')
  # print('number of edges', ls.shape)
  cubes = get_cubes_from_graph(G_randadded, Grid_Number)
  sio.savemat(full_path, {'cubes': cubes*cube_size/(Grid_Number), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})
  pass


def Periodic_Cube_Run(p1, p2, bx, by, bz, Grid_Number, cube_size, NameData):
  "generate periodic cube based unit cell"

  # Set the output file name
  filename = NameData+".mat"
  # full_path = path + "\\" + filename
  full_path = os.path.join(path, filename)

  # G_complete = get_complete_graph_cube(n_grid=Grid_Number)
  G_interior,G_exterior_xy_0,G_exterior_yz_0,G_exterior_xz_0 = get_complete_graph_cube_with_Boundaries(Grid_Number)

  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=10:
      G_interior_randsub = get_random_subgraph(G_interior,p=p1)
      print("Interior")
      G_interior_span_tree = get_random_spanning_tree(G_interior_randsub)
      if G_interior_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1
  G_interior_randadded = add_random_complement_graph(G_interior_span_tree,G_interior_randsub,p=p2)

  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=10:
      G_exterior_xy_0_randsub = get_random_subgraph(G_exterior_xy_0,p=p1)
      print("G_exterior_xy_0")
      G_exterior_xy_0_span_tree = get_random_spanning_tree(G_exterior_xy_0_randsub)
      if G_exterior_xy_0_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1
  G_exterior_xy_0_randadded = add_random_complement_graph(G_exterior_xy_0_span_tree,G_exterior_xy_0_randsub,p=p2)

  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=10:
      G_exterior_yz_0_randsub = get_random_subgraph(G_exterior_yz_0,p=p1)
      print("G_exterior_yz_0")
      G_exterior_yz_0_span_tree = get_random_spanning_tree(G_exterior_yz_0_randsub)
      if G_exterior_yz_0_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1
  G_exterior_yz_0_randadded = add_random_complement_graph(G_exterior_yz_0_span_tree,G_exterior_yz_0_randsub,p=p2)

  SpanningTreeSuccessful = True
  NumberOfIteration = 1
  while SpanningTreeSuccessful and NumberOfIteration<=10:
      G_exterior_xz_0_randsub = get_random_subgraph(G_exterior_xz_0,p=p1)
      print("G_exterior_xz_0")
      G_exterior_xz_0_span_tree = get_random_spanning_tree(G_exterior_xz_0_randsub)
      if G_exterior_xz_0_span_tree!= None:
          SpanningTreeSuccessful = False
      NumberOfIteration+=1
  G_exterior_xz_0_randadded = add_random_complement_graph(G_exterior_xz_0_span_tree,G_exterior_xz_0_randsub,p=p2)

  cubes = get_cubes_from_graph_with_boundaries(Grid_Number, G_interior_randadded, G_exterior_xy_0_randadded, G_exterior_yz_0_randadded, G_exterior_xz_0_randadded)

  sio.savemat(full_path, {'cubes': cubes*cube_size/(Grid_Number), 'p1': p1, 'p2': p2, 'bx': bx, 'by': by, 'bz': bz, 'Grid_Number': Grid_Number})

  pass

##----------------------------------Calling the codes--------------------------------------##

# Set the path to the directory where you want to save the file
path = r"C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial"
# path = r"/home/rkm41/GraphDataZ"

parser = argparse.ArgumentParser()
parser.add_argument('--NameData')
parser.add_argument('--p1', type=float)
parser.add_argument('--p2', type=float)
parser.add_argument('--bx', type=float)
parser.add_argument('--by', type=float)
parser.add_argument('--bz', type=float)
parser.add_argument('--Grid_Number', type=int)
parser.add_argument('--cube_size', type=int)
parser.add_argument('--Base_Element', choices=['Rod', 'Triangle', 'Pyramid', 'Cube', 'PeriodicCube'])

args = parser.parse_args()

print(args.p1)
print(args.p2)
print(args.bx)
print(args.by)
print(args.bz)
print(args.Grid_Number)
print(args.cube_size)
print(args.Base_Element)

start_time = time.time()
if args.Base_Element == 'Rod':
  Rod_Run(args.p1, args.p2, args.bx, args.by, args.bz, args.Grid_Number, args.cube_size, args.NameData)
elif args.Base_Element == 'Triangle':
  Triangle_Run(args.p1, args.p2, args.bx, args.by, args.bz, args.Grid_Number, args.cube_size, args.NameData)
elif args.Base_Element == 'Pyramid':
  Pyramid_Run(args.p1, args.p2, args.bx, args.by, args.bz, args.Grid_Number, args.cube_size, args.NameData)
elif args.Base_Element == 'Cube':
  Cube_Run(args.p1, args.p2, args.bx, args.by, args.bz, args.Grid_Number, args.cube_size, args.NameData)
elif args.Base_Element == 'PeriodicCube':
  Periodic_Cube_Run(args.p1, args.p2, args.bx, args.by, args.bz, args.Grid_Number, args.cube_size, args.NameData)
# else:
#   error_msg = f"Invalid value '{args.Base_Element}' for --Base_Element. Expected 'Rod' or 'Triangle' or 'Pyramid' or 'Cube' or 'PeriodicCube'."
#   raise ValueError(error_msg)

# Stop the timer
end_time = time.time()

# Calculate the duration
runtime = end_time - start_time

print(f"Runtime: {runtime} seconds")
