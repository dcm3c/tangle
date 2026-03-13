#!/usr/bin/env python3
"""
showmesh.py - Visualize Triangle-format meshes (.node/.ele/.poly files)

Usage:
    python showmesh.py base              # show base.1.node + base.1.ele
    python showmesh.py base.node         # show base.node + base.ele
    python showmesh.py base.poly         # show input PSLG (segments + holes)
    python showmesh.py -p base           # show input base.poly overlaid on mesh
    python showmesh.py -n base           # label vertex numbers
    python showmesh.py -e base           # label element numbers
    python showmesh.py -r base           # color by region attribute
    python showmesh.py -o out.png base   # save to file instead of showing
"""

import sys
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
import numpy as np


def read_node(fn):
    pts = {}
    with open(fn) as f:
        parts = f.readline().split()
        n = int(parts[0])
        n_attr = int(parts[1]) if len(parts) > 1 else 0
        n_mark = int(parts[2]) if len(parts) > 2 else 0
        has_mark = int(parts[3]) if len(parts) > 3 else 0
        for _ in range(n):
            parts = f.readline().split()
            if not parts:
                continue
            idx = int(parts[0])
            x, y = float(parts[1]), float(parts[2])
            marker = int(parts[3 + n_attr]) if has_mark and len(parts) > 3 + n_attr else 0
            pts[idx] = (x, y, marker)
    return pts


def read_ele(fn):
    tris = []
    attrs = []
    with open(fn) as f:
        parts = f.readline().split()
        n = int(parts[0])
        n_per = int(parts[1]) if len(parts) > 1 else 3
        n_attr = int(parts[2]) if len(parts) > 2 else 0
        for _ in range(n):
            parts = f.readline().split()
            if not parts:
                continue
            v = [int(parts[j]) for j in range(1, n_per + 1)]
            attr = int(float(parts[n_per + 1])) if n_attr > 0 and len(parts) > n_per + 1 else 0
            tris.append(v)
            attrs.append(attr)
    return tris, attrs


def read_poly(fn):
    """Read a .poly file, returning vertices, segments, and holes."""
    verts = {}
    segs = []
    holes = []
    with open(fn) as f:
        # Vertices
        parts = f.readline().split()
        nv = int(parts[0])
        ndim = int(parts[1]) if len(parts) > 1 else 2
        n_attr = int(parts[2]) if len(parts) > 2 else 0
        has_mark = int(parts[3]) if len(parts) > 3 else 0
        for _ in range(nv):
            parts = f.readline().split()
            if not parts:
                continue
            idx = int(parts[0])
            x, y = float(parts[1]), float(parts[2])
            verts[idx] = (x, y)
        # Segments
        parts = f.readline().split()
        ns = int(parts[0])
        seg_has_mark = int(parts[1]) if len(parts) > 1 else 0
        for _ in range(ns):
            parts = f.readline().split()
            if not parts:
                continue
            v0, v1 = int(parts[1]), int(parts[2])
            segs.append((v0, v1))
        # Holes
        parts = f.readline().split()
        if parts:
            nh = int(parts[0])
            for _ in range(nh):
                parts = f.readline().split()
                if not parts:
                    continue
                hx, hy = float(parts[1]), float(parts[2])
                holes.append((hx, hy))
    return verts, segs, holes


def find_files(base):
    """Given a base name, find the .node and .ele files."""
    # Strip known extensions
    for ext in ('.node', '.ele', '.poly', '.edge', '.neigh'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    # Try base.1.node first (output), then base.node (input)
    for suffix in ['.1', '']:
        node = base + suffix + '.node'
        ele = base + suffix + '.ele'
        if os.path.exists(node) and os.path.exists(ele):
            return base, node, ele

    return base, None, None


def main():
    parser = argparse.ArgumentParser(description='Visualize Triangle-format meshes')
    parser.add_argument('input', help='Base name, .node, .ele, or .poly file')
    parser.add_argument('-p', '--poly', action='store_true',
                        help='Overlay input .poly segments')
    parser.add_argument('-n', '--node-labels', action='store_true',
                        help='Label vertex numbers')
    parser.add_argument('-e', '--ele-labels', action='store_true',
                        help='Label element numbers')
    parser.add_argument('-r', '--regions', action='store_true',
                        help='Color triangles by region attribute')
    parser.add_argument('-o', '--output', metavar='FILE',
                        help='Save to file instead of showing')
    parser.add_argument('--no-mesh', action='store_true',
                        help='Hide mesh edges (useful with -p to show only PSLG)')
    args = parser.parse_args()

    inp = args.input

    # If input is a .poly file and no mesh exists, just show the PSLG
    poly_only = False
    if inp.endswith('.poly') and not os.path.exists(inp.replace('.poly', '.ele')):
        base_name = inp[:-5]
        _, node_file, ele_file = find_files(base_name)
        if node_file is None:
            poly_only = True

    if poly_only:
        verts, segs, holes = read_poly(inp)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        ax.set_aspect('equal')
        ax.set_title(os.path.basename(inp))

        # Draw segments
        lines = []
        for v0, v1 in segs:
            if v0 in verts and v1 in verts:
                lines.append([verts[v0][:2], verts[v1][:2]])
        if lines:
            lc = LineCollection(lines, colors='blue', linewidths=1.5)
            ax.add_collection(lc)

        # Draw vertices
        xs = [v[0] for v in verts.values()]
        ys = [v[1] for v in verts.values()]
        ax.plot(xs, ys, 'k.', markersize=3)

        # Draw holes
        for hx, hy in holes:
            ax.plot(hx, hy, 'rx', markersize=8, markeredgewidth=2)

        if args.node_labels:
            for idx, (x, y) in verts.items():
                ax.annotate(str(idx), (x, y), fontsize=6, ha='center', va='bottom')

        ax.autoscale()
        ax.margins(0.05)

        if args.output:
            fig.savefig(args.output, dpi=150, bbox_inches='tight')
            print(f"Saved to {args.output}")
        else:
            plt.show()
        return

    # Find mesh files
    base_name, node_file, ele_file = find_files(inp)
    if node_file is None:
        print(f"Error: cannot find .node/.ele files for '{inp}'", file=sys.stderr)
        sys.exit(1)

    pts = read_node(node_file)
    tris, attrs = read_ele(ele_file)

    # Build index mapping (handles both 0-based and 1-based)
    min_idx = min(pts.keys())
    idx_list = sorted(pts.keys())
    idx_to_pos = {idx: i for i, idx in enumerate(idx_list)}
    x = np.array([pts[i][0] for i in idx_list])
    y = np.array([pts[i][1] for i in idx_list])

    # Remap triangle indices
    tri_array = np.array([[idx_to_pos[v] for v in t] for t in tris])

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.set_aspect('equal')
    ax.set_title(f"{os.path.basename(node_file)} ({len(tris)} triangles, {len(pts)} vertices)")

    if args.regions and any(a != 0 for a in attrs):
        # Color by region
        triang = mtri.Triangulation(x, y, tri_array)
        ax.tripcolor(triang, facecolors=np.array(attrs, dtype=float),
                     cmap='Set3', edgecolors='gray', linewidth=0.3)
    elif not args.no_mesh:
        # Plain mesh
        ax.triplot(x, y, tri_array, 'k-', linewidth=0.3)

    if args.node_labels:
        for idx in idx_list:
            px, py, _ = pts[idx]
            ax.annotate(str(idx), (px, py), fontsize=5, ha='center', va='bottom',
                        color='blue')

    if args.ele_labels:
        for i, t in enumerate(tris):
            cx = np.mean([pts[v][0] for v in t])
            cy = np.mean([pts[v][1] for v in t])
            label = i + min(int(p[0]) for p in [open(ele_file).readline().split()])
            ax.annotate(str(i + min_idx), (cx, cy), fontsize=4, ha='center',
                        va='center', color='red')

    # Overlay .poly segments if requested
    if args.poly:
        poly_file = base_name + '.poly'
        if not os.path.exists(poly_file):
            print(f"Warning: {poly_file} not found", file=sys.stderr)
        else:
            verts, segs, holes = read_poly(poly_file)
            lines = []
            for v0, v1 in segs:
                if v0 in verts and v1 in verts:
                    lines.append([verts[v0][:2], verts[v1][:2]])
            if lines:
                lc = LineCollection(lines, colors='blue', linewidths=1.5, zorder=5)
                ax.add_collection(lc)
            for hx, hy in holes:
                ax.plot(hx, hy, 'rx', markersize=8, markeredgewidth=2, zorder=5)

    ax.margins(0.02)

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
