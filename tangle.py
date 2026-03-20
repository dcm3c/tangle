#!/usr/bin/env python3
"""
tangle
A 2D Delaunay triangulation tool compatible with Shewchuk's Triangle format.

Author: David Meeker
Generated with the assistance of Claude Code Opus 4.6
Translated from tangle.cpp

Version 0.2
17 Mar 2026

Supports: -p -P -j -q -e -A -a -z -Q -I -Y options

Usage: python tangle.py [switches] inputfile

File format compatibility:
  .node  - vertex files
  .poly  - PSLG (Planar Straight-Line Graph) files
  .ele   - triangle output files
  .edge  - edge output files
  .neigh - neighbor output files

MIT License

Copyright (c) 2026 David Meeker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import math
import time
import bisect
from collections import deque
from collections import deque
from dataclasses import dataclass, field
from decimal import Decimal
from heapq import heappush, heappop

# Base relative tolerance.  EPS is scaled to the bounding-box diagonal in
# build_delaunay(); EPS_SCALE controls that ratio and is also used directly
# for dimensionless comparisons.
EPS_SCALE = 1e-10
EPS = EPS_SCALE


# ============================================================
# Options
# ============================================================

@dataclass
class Options:
    pslg: bool = False
    no_poly_out: bool = False
    jettison: bool = False
    quality: bool = False
    min_angle: float = 20.0
    edges: bool = False
    regions: bool = False
    area_limit: bool = False
    max_area: float = -1.0
    zero_indexed: bool = False
    quiet: bool = False
    suppress_iter: bool = False
    no_holes: bool = False
    no_steiner: bool = False
    neighbors: bool = False
    clean_pslg: bool = False
    clean_tol: float = -1.0  # -1 = auto (bboxDiag * 1e-6)
    first_index: int = 1


# ============================================================
# Geometry primitives
# ============================================================

@dataclass
class Point:
    x: float
    y: float
    id: int = 0
    marker: int = 0
    attribs: list = field(default_factory=list)
    lfs: float = -1.0  # local feature size constraint (-1 = none)


@dataclass
class Segment:
    v0: int
    v1: int
    marker: int = 0
    lfs: float = -1.0  # per-segment LFS constraint (-1 = none)
    pbc_partner: int = -1  # index of paired PBC segment (-1 = none)
    pbc_type: int = -1     # 0=periodic, 1=anti-periodic (-1 = none)


@dataclass
class Hole:
    x: float
    y: float


@dataclass
class Region:
    x: float
    y: float
    attrib: float = 0.0
    max_area: float = -1.0


@dataclass
class Triangle:
    v: list = field(default_factory=lambda: [-1, -1, -1])
    neighbors: list = field(default_factory=lambda: [-1, -1, -1])
    region_attrib: float = 0.0
    region_max_area: float = -1.0
    marker: int = 0
    generation: int = 0


def orient2d_xy(ax, ay, bx, by, cx, cy):
    """Lightweight orient2d for triangle walks — no exact fallback needed because
    the walk is intrinsically robust to occasional sign errors (it self-corrects
    or falls through to the linear scan).  Also avoids constructing a Point."""
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def orient2d(a, b, c):
    Ax, Ay = b.x - a.x, b.y - a.y
    Bx, By = c.x - a.x, c.y - a.y
    det = Ax * By - Bx * Ay

    # Error bound: det terms are O(M²), check if result is well above fp noise
    M = max(abs(Ax), abs(Ay), abs(Bx), abs(By))
    if abs(det) > 6.662e-16 * M * M:
        return det

    # Fall back to exact arithmetic via Decimal
    fAx = Decimal(b.x) - Decimal(a.x)
    fAy = Decimal(b.y) - Decimal(a.y)
    fBx = Decimal(c.x) - Decimal(a.x)
    fBy = Decimal(c.y) - Decimal(a.y)
    return float(fAx * fBy - fBx * fAy)


def in_circle(a, b, c, d):
    # 3x3 formulation: subtract d to eliminate one row
    Ax, Ay = a.x - d.x, a.y - d.y
    Bx, By = b.x - d.x, b.y - d.y
    Cx, Cy = c.x - d.x, c.y - d.y
    As = Ax * Ax + Ay * Ay
    Bs = Bx * Bx + By * By
    Cs = Cx * Cx + Cy * Cy

    det = (Ax * (By * Cs - Cy * Bs)
         - Ay * (Bx * Cs - Cx * Bs)
         + As * (Bx * Cy - Cx * By))

    # Error bound: det terms are O(M⁴), check if result is well above fp noise
    M = max(abs(Ax), abs(Ay), abs(Bx), abs(By),
            abs(Cx), abs(Cy))
    M4 = M * M * M * M
    if abs(det) > 1.333e-14 * M4:
        return det

    # Fall back to exact arithmetic via Decimal (standard 3x3 formulation)
    ax = Decimal(a.x) - Decimal(d.x)
    ay = Decimal(a.y) - Decimal(d.y)
    bx = Decimal(b.x) - Decimal(d.x)
    by = Decimal(b.y) - Decimal(d.y)
    cx = Decimal(c.x) - Decimal(d.x)
    cy = Decimal(c.y) - Decimal(d.y)
    cs = cx * cx + cy * cy
    bs = bx * bx + by * by
    fdet = (ax * (by * cs - cy * bs)
          - ay * (bx * cs - cx * bs)
          + (ax * ax + ay * ay) * (bx * cy - cx * by))
    return float(fdet)


def tri_area(a, b, c):
    return 0.5 * abs(orient2d(a, b, c))


def vertex_angle(a, b, c):
    abx, aby = b.x - a.x, b.y - a.y
    acx, acy = c.x - a.x, c.y - a.y
    dot = abx * acx + aby * acy
    cross = abx * acy - aby * acx
    return math.atan2(abs(cross), dot) * 180.0 / math.pi


def min_angle(a, b, c):
    la = math.sqrt((b.x-c.x)*(b.x-c.x)+(b.y-c.y)*(b.y-c.y))
    lb = math.sqrt((a.x-c.x)*(a.x-c.x)+(a.y-c.y)*(a.y-c.y))
    lc = math.sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))
    if la < EPS or lb < EPS or lc < EPS:
        return 0.0
    cos_a = max(-1.0, min(1.0, (lb * lb + lc * lc - la * la) / (2 * lb * lc)))
    cos_b = max(-1.0, min(1.0, (la * la + lc * lc - lb * lb) / (2 * la * lc)))
    cos_c = max(-1.0, min(1.0, (la * la + lb * lb - lc * lc) / (2 * la * lb)))
    return min(math.acos(cos_a), math.acos(cos_b), math.acos(cos_c)) * 180.0 / math.pi


def circumcenter(a, b, c):
    ax, ay = b.x - a.x, b.y - a.y
    bx, by = c.x - a.x, c.y - a.y
    D = 2.0 * (ax * by - ay * bx)
    if abs(D) < EPS:
        return (a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0
    ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / D
    uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / D
    return a.x + ux, a.y + uy


def point_on_segment(p, a, b):
    if orient2d(a, b, p) != 0:
        return False
    dx, dy = b.x - a.x, b.y - a.y
    if abs(dx) > EPS:
        t = (p.x - a.x) / dx
    elif abs(dy) > EPS:
        t = (p.y - a.y) / dy
    else:
        return False
    return -EPS_SCALE <= t <= 1 + EPS_SCALE  # dimensionless: fixed tolerance


def edge_key(a, b):
    return (min(a, b), max(a, b))


# ============================================================
# Mesh structure
# ============================================================

@dataclass
class PBCNodePair:
    node_a: int
    node_b: int
    type: int  # 0=periodic, 1=anti-periodic


class Mesh:
    def __init__(self):
        self.vertices = []
        self.triangles = []
        self.segments = []
        self.holes = []
        self.regions = []
        self.edges = []
        self.pbc_pairs = []      # list of PBCNodePair, filled after refinement
        self.pbc_twin = {}       # dict: node -> twin node (bidirectional)
        self.pbc_node_type = {}  # dict: node -> pbc type (0 or 1)

    def rebuild_adjacency(self):
        nt = len(self.triangles)
        edge_to_tri_edge = {}
        for t in self.triangles:
            t.neighbors = [-1, -1, -1]
        for i in range(nt):
            t = self.triangles[i]
            if t.v[0] < 0:
                continue
            for j in range(3):
                a, b = t.v[j], t.v[(j + 1) % 3]
                mn, mx = min(a, b), max(a, b)
                key = (mn << 32) | mx
                if key in edge_to_tri_edge:
                    packed = edge_to_tri_edge[key]
                    ot, oe = packed >> 2, packed & 3
                    self.triangles[i].neighbors[j] = ot
                    self.triangles[ot].neighbors[oe] = i
                    del edge_to_tri_edge[key]
                else:
                    edge_to_tri_edge[key] = (i << 2) | j

    def locate_triangle(self, px, py, hint=0):
        tris = self.triangles
        if not tris:
            return -1
        n = len(tris)
        if hint < 0 or hint >= n:
            hint = 0
        cur = hint
        verts = self.vertices
        for _ in range(n * 2):
            t = tris[cur]
            if t.v[0] < 0:
                cur = (cur + 1) % n
                continue
            v0, v1, v2 = verts[t.v[0]], verts[t.v[1]], verts[t.v[2]]
            o0 = orient2d_xy(v0.x, v0.y, v1.x, v1.y, px, py)
            o1 = orient2d_xy(v1.x, v1.y, v2.x, v2.y, px, py)
            o2 = orient2d_xy(v2.x, v2.y, v0.x, v0.y, px, py)
            if o0 >= -EPS and o1 >= -EPS and o2 >= -EPS:
                return cur
            worst, wv = -1, 0
            if o0 < wv:
                wv, worst = o0, 0
            if o1 < wv:
                wv, worst = o1, 1
            if o2 < wv:
                wv, worst = o2, 2
            if worst >= 0 and t.neighbors[worst] >= 0:
                cur = t.neighbors[worst]
            else:
                break
        # Fallback: linear scan
        for i in range(n):
            t = tris[i]
            if t.v[0] < 0:
                continue
            v0, v1, v2 = verts[t.v[0]], verts[t.v[1]], verts[t.v[2]]
            if (orient2d_xy(v0.x, v0.y, v1.x, v1.y, px, py) >= -EPS and
                orient2d_xy(v1.x, v1.y, v2.x, v2.y, px, py) >= -EPS and
                orient2d_xy(v2.x, v2.y, v0.x, v0.y, px, py) >= -EPS):
                return i
        return -1

    def has_edge(self, u, v, v2t):
        for ti in v2t[u]:
            t = self.triangles[ti]
            if t.v[0] < 0:
                continue
            for j in range(3):
                if (t.v[j] == u and t.v[(j + 1) % 3] == v) or \
                   (t.v[j] == v and t.v[(j + 1) % 3] == u):
                    return True
        return False


# ============================================================
# PSLG cleanup (-C flag)
# ============================================================

def clean_pslg(mesh, tol, quiet):
    pts = mesh.vertices
    segs = mesh.segments
    nv = len(pts)

    # Compute auto tolerance if not specified
    if tol < 0:
        min_x = min(p.x for p in pts)
        max_x = max(p.x for p in pts)
        min_y = min(p.y for p in pts)
        max_y = max(p.y for p in pts)
        tol = math.sqrt((max_x - min_x)**2 + (max_y - min_y)**2) * 1e-6
    tol2 = tol * tol

    merged_nodes = 0
    split_segs = 0
    removed_segs = 0

    # 1. Merge near-duplicate nodes
    remap = list(range(nv))
    for i in range(nv):
        if remap[i] != i:
            continue
        for j in range(i + 1, nv):
            if remap[j] != j:
                continue
            dx = pts[j].x - pts[i].x
            dy = pts[j].y - pts[i].y
            if dx * dx + dy * dy < tol2:
                remap[j] = i
                merged_nodes += 1

    # Apply remap to segments
    for s in segs:
        s.v0 = remap[s.v0]
        s.v1 = remap[s.v1]

    # 2. Remove degenerate segments
    good = [s for s in segs if s.v0 != s.v1]
    removed_segs += len(segs) - len(good)
    mesh.segments = good
    segs = mesh.segments

    # 3. Remove duplicate segments
    seen = set()
    good = []
    for s in segs:
        k = edge_key(s.v0, s.v1)
        if k not in seen:
            seen.add(k)
            good.append(s)
        else:
            removed_segs += 1
    mesh.segments = good
    segs = mesh.segments

    # 4. Split segments at near-coincident nodes
    changed = True
    while changed:
        changed = False
        for si in range(len(segs)):
            v0, v1 = segs[si].v0, segs[si].v1
            ax, ay = pts[v0].x, pts[v0].y
            bx, by = pts[v1].x, pts[v1].y
            dx, dy = bx - ax, by - ay
            len2 = dx * dx + dy * dy
            if len2 < tol2:
                continue
            sqrt_len = math.sqrt(len2)
            for i in range(len(pts)):
                if i == v0 or i == v1:
                    continue
                if i < nv and remap[i] != i:
                    continue
                px, py = pts[i].x - ax, pts[i].y - ay
                t = (dx * px + dy * py) / len2
                if t <= tol / sqrt_len or t >= 1.0 - tol / sqrt_len:
                    continue
                cross = dx * py - dy * px
                if cross * cross > tol2 * len2:
                    continue
                # Point i is on segment si — split it
                s2 = Segment(i, segs[si].v1, segs[si].marker)
                segs[si] = Segment(segs[si].v0, i, segs[si].marker)
                segs.append(s2)
                split_segs += 1
                changed = True
                break

    # 5. Split segments at segment-segment intersections
    changed = True
    while changed:
        changed = False
        for i in range(len(segs)):
            if changed:
                break
            ax, ay = pts[segs[i].v0].x, pts[segs[i].v0].y
            bx, by = pts[segs[i].v1].x, pts[segs[i].v1].y
            for j in range(i + 1, len(segs)):
                if changed:
                    break
                # Skip if they share an endpoint
                if (segs[i].v0 in (segs[j].v0, segs[j].v1) or
                    segs[i].v1 in (segs[j].v0, segs[j].v1)):
                    continue
                cx, cy = pts[segs[j].v0].x, pts[segs[j].v0].y
                d2x, d2y = pts[segs[j].v1].x, pts[segs[j].v1].y
                d1x, d1y = bx - ax, by - ay
                ddx, ddy = d2x - cx, d2y - cy
                denom = d1x * ddy - d1y * ddx
                if abs(denom) < 1e-15:
                    continue
                t = ((cx - ax) * ddy - (cy - ay) * ddx) / denom
                u = ((cx - ax) * d1y - (cy - ay) * d1x) / denom
                if t <= 0 or t >= 1 or u <= 0 or u >= 1:
                    continue
                ix, iy = ax + t * d1x, ay + t * d1y
                # Check if near existing node
                near_node = -1
                for k in range(len(pts)):
                    ex, ey = pts[k].x - ix, pts[k].y - iy
                    if ex * ex + ey * ey < tol2:
                        near_node = k
                        break
                if near_node < 0:
                    near_node = len(pts)
                    pts.append(Point(ix, iy, near_node, 0))
                # Split both segments
                s_i2 = Segment(near_node, segs[i].v1, segs[i].marker)
                segs[i] = Segment(segs[i].v0, near_node, segs[i].marker)
                segs.append(s_i2)
                s_j2 = Segment(near_node, segs[j].v1, segs[j].marker)
                segs[j] = Segment(segs[j].v0, near_node, segs[j].marker)
                segs.append(s_j2)
                split_segs += 2
                changed = True

    # Compact: remove merged-away vertices and re-index
    if merged_nodes > 0:
        new_idx = [-1] * len(pts)
        new_pts = []
        for i in range(len(pts)):
            target = remap[i] if i < nv else i
            if target != i:
                continue
            new_idx[i] = len(new_pts)
            pts[i].id = len(new_pts)
            new_pts.append(pts[i])
        for i in range(nv):
            if remap[i] != i:
                new_idx[i] = new_idx[remap[i]]
        for i in range(nv, len(pts)):
            if new_idx[i] < 0:
                new_idx[i] = len(new_pts)
                pts[i].id = len(new_pts)
                new_pts.append(pts[i])
        for s in segs:
            s.v0 = new_idx[s.v0]
            s.v1 = new_idx[s.v1]
        mesh.vertices = new_pts

    if not quiet and (merged_nodes > 0 or split_segs > 0 or removed_segs > 0):
        print(f"  PSLG cleanup: merged {merged_nodes} nodes, split "
              f"{split_segs} segments, removed {removed_segs} duplicates",
              file=sys.stderr)


# ============================================================
# Bowyer-Watson Delaunay triangulation
# ============================================================

def build_delaunay(mesh):
    pts = mesh.vertices
    n = len(pts)

    min_x = min(p.x for p in pts)
    max_x = max(p.x for p in pts)
    min_y = min(p.y for p in pts)
    max_y = max(p.y for p in pts)
    dx, dy = max_x - min_x, max_y - min_y
    d = max(dx, dy) * 10.0

    # Scale EPS to the bounding-box diagonal so tolerances stay meaningful
    # regardless of coordinate magnitude.
    global EPS
    bbox_diag = math.sqrt(dx*dx + dy*dy)
    EPS = bbox_diag * EPS_SCALE
    pts.append(Point(min_x - d, min_y - 3 * d, n, 0))
    pts.append(Point(min_x + 3 * d, min_y - d, n + 1, 0))
    pts.append(Point(min_x - d, min_y + 3 * d, n + 2, 0))

    tris = mesh.triangles
    tris.clear()
    tris.append(Triangle(v=[n, n + 1, n + 2], neighbors=[-1, -1, -1]))

    order = sorted(range(n), key=lambda i: (pts[i].x, pts[i].y))

    last_tri = 0
    visit_epoch = []
    epoch = 0

    for ii in range(n):
        pidx = order[ii]
        p = pts[pidx]

        # Locate containing triangle by walking
        start_tri = -1
        cur = min(last_tri, len(tris) - 1)
        for _ in range(len(tris) * 2):
            t = tris[cur]
            if t.v[0] < 0:
                cur = 0
                continue
            o0 = orient2d(pts[t.v[0]], pts[t.v[1]], p)
            o1 = orient2d(pts[t.v[1]], pts[t.v[2]], p)
            o2 = orient2d(pts[t.v[2]], pts[t.v[0]], p)
            if o0 >= -EPS and o1 >= -EPS and o2 >= -EPS:
                start_tri = cur
                break
            w, wv = -1, 0
            if o0 < wv:
                wv, w = o0, 0
            if o1 < wv:
                wv, w = o1, 1
            if o2 < wv:
                wv, w = o2, 2
            if w >= 0 and t.neighbors[w] >= 0:
                cur = t.neighbors[w]
            else:
                break
        if start_tri < 0:
            for i in range(len(tris)):
                t = tris[i]
                if t.v[0] < 0:
                    continue
                if (orient2d(pts[t.v[0]], pts[t.v[1]], p) >= -EPS and
                    orient2d(pts[t.v[1]], pts[t.v[2]], p) >= -EPS and
                    orient2d(pts[t.v[2]], pts[t.v[0]], p) >= -EPS):
                    start_tri = i
                    break
        if start_tri < 0:
            continue

        # BFS cavity (epoch-based visited to avoid per-point allocation)
        epoch += 1
        while len(visit_epoch) < len(tris):
            visit_epoch.append(0)
        bad = []
        stk = [start_tri]
        while stk:
            cur = stk.pop()
            if cur < 0 or cur >= len(tris) or visit_epoch[cur] == epoch:
                continue
            visit_epoch[cur] = epoch
            t = tris[cur]
            if t.v[0] < 0:
                continue
            if in_circle(pts[t.v[0]], pts[t.v[1]], pts[t.v[2]], p) > 0:
                bad.append(cur)
                for j in range(3):
                    if t.neighbors[j] >= 0 and visit_epoch[t.neighbors[j]] != epoch:
                        stk.append(t.neighbors[j])

        if not bad:
            bad.append(start_tri)

        # Collect boundary polygon
        bad_set = set(bad)
        poly = []  # list of (v0, v1, outer_tri, outer_local_edge)
        for bi in bad:
            t = tris[bi]
            for j in range(3):
                nb = t.neighbors[j]
                if nb < 0 or nb not in bad_set:
                    oe = -1
                    if nb >= 0:
                        for k in range(3):
                            if tris[nb].neighbors[k] == bi:
                                oe = k
                                break
                    poly.append((t.v[j], t.v[(j + 1) % 3], nb, oe))

        n_new = len(poly)
        bad.sort()
        slots = []
        for i in range(min(n_new, len(bad))):
            slots.append(bad[i])
        while len(slots) < n_new:
            slots.append(len(tris))
            tris.append(Triangle())
        for i in range(n_new, len(bad)):
            tris[bad[i]].v = [-1, -1, -1]
            tris[bad[i]].neighbors = [-1, -1, -1]

        for i in range(n_new):
            s = slots[i]
            tris[s].v = [poly[i][0], poly[i][1], pidx]
            tris[s].neighbors = [-1, -1, -1]
            tris[s].region_attrib = 0
            tris[s].region_max_area = -1
            tris[s].marker = 0
            if orient2d(pts[tris[s].v[0]], pts[tris[s].v[1]], pts[tris[s].v[2]]) < 0:
                tris[s].v[0], tris[s].v[1] = tris[s].v[1], tris[s].v[0]

        # Wire outer adjacency
        for i in range(n_new):
            s = slots[i]
            if poly[i][2] < 0:
                continue
            for j in range(3):
                a, b = tris[s].v[j], tris[s].v[(j + 1) % 3]
                if (a == poly[i][0] and b == poly[i][1]) or \
                   (a == poly[i][1] and b == poly[i][0]):
                    tris[s].neighbors[j] = poly[i][2]
                    if poly[i][3] >= 0:
                        tris[poly[i][2]].neighbors[poly[i][3]] = s
                    break

        # Wire internal adjacency
        pidx_edges = {}
        for i in range(n_new):
            s = slots[i]
            for j in range(3):
                a, b = tris[s].v[j], tris[s].v[(j + 1) % 3]
                if a != pidx and b != pidx:
                    continue
                other = b if a == pidx else a
                if other in pidx_edges:
                    os, oe = pidx_edges[other]
                    tris[s].neighbors[j] = os
                    tris[os].neighbors[oe] = s
                    del pidx_edges[other]
                else:
                    pidx_edges[other] = (s, j)
        last_tri = slots[0]

    # Remove dead triangles and super-triangle vertices
    keep = [False] * len(tris)
    for i in range(len(tris)):
        if tris[i].v[0] < 0:
            continue
        if tris[i].v[0] >= n or tris[i].v[1] >= n or tris[i].v[2] >= n:
            continue
        keep[i] = True
    remap = [-1] * len(tris)
    new_tris = []
    for i in range(len(tris)):
        if keep[i]:
            remap[i] = len(new_tris)
            new_tris.append(tris[i])
    for t in new_tris:
        for j in range(3):
            t.neighbors[j] = remap[t.neighbors[j]] if t.neighbors[j] >= 0 else -1
    mesh.triangles = new_tris
    del pts[n:]

    # Detect and re-insert orphaned vertices
    mesh.rebuild_adjacency()
    tris = mesh.triangles
    orphan_visit_epoch = []
    orphan_epoch = 0
    for retry in range(10):
        in_mesh = [False] * n
        for t in tris:
            if t.v[0] >= 0:
                for j in range(3):
                    if t.v[j] < n:
                        in_mesh[t.v[j]] = True
        orphans = [i for i in range(n) if not in_mesh[i]]
        if not orphans:
            break

        for pidx in orphans:
            p = pts[pidx]
            contain_tri = mesh.locate_triangle(p.x, p.y)
            if contain_tri < 0:
                continue

            # Skip if nearly coincident
            too_close = False
            for j in range(3):
                dx2 = pts[tris[contain_tri].v[j]].x - p.x
                dy2 = pts[tris[contain_tri].v[j]].y - p.y
                if dx2 * dx2 + dy2 * dy2 < EPS * EPS:
                    too_close = True
                    break
            if too_close:
                continue

            # Bowyer-Watson insertion with strict inCircle to avoid over-expanding
            # the cavity and orphaning nearby vertices
            orphan_epoch += 1
            while len(orphan_visit_epoch) < len(tris):
                orphan_visit_epoch.append(0)
            bad = []
            stk = [contain_tri]
            while stk:
                cur = stk.pop()
                if cur < 0 or cur >= len(tris) or orphan_visit_epoch[cur] == orphan_epoch:
                    continue
                orphan_visit_epoch[cur] = orphan_epoch
                t = tris[cur]
                if t.v[0] < 0:
                    continue
                if in_circle(pts[t.v[0]], pts[t.v[1]], pts[t.v[2]], p) > 0:
                    bad.append(cur)
                    for j in range(3):
                        if t.neighbors[j] >= 0 and orphan_visit_epoch[t.neighbors[j]] != orphan_epoch:
                            stk.append(t.neighbors[j])
            if not bad:
                bad.append(contain_tri)

            bad_set = set(bad)
            poly = []
            for bi in bad:
                t = tris[bi]
                for j in range(3):
                    nb = t.neighbors[j]
                    if nb < 0 or nb not in bad_set:
                        oe = -1
                        if nb >= 0:
                            for k in range(3):
                                if tris[nb].neighbors[k] == bi:
                                    oe = k
                                    break
                        poly.append((t.v[j], t.v[(j + 1) % 3], nb, oe))

            n_new = len(poly)
            bad.sort()
            slots = []
            for i in range(min(n_new, len(bad))):
                slots.append(bad[i])
            while len(slots) < n_new:
                slots.append(len(tris))
                tris.append(Triangle())
            for i in range(n_new, len(bad)):
                tris[bad[i]].v = [-1, -1, -1]
                tris[bad[i]].neighbors = [-1, -1, -1]

            for i in range(n_new):
                s = slots[i]
                tris[s].v = [poly[i][0], poly[i][1], pidx]
                tris[s].neighbors = [-1, -1, -1]
                tris[s].region_attrib = 0
                tris[s].region_max_area = -1
                tris[s].marker = 0
                if orient2d(pts[tris[s].v[0]], pts[tris[s].v[1]], pts[tris[s].v[2]]) < 0:
                    tris[s].v[0], tris[s].v[1] = tris[s].v[1], tris[s].v[0]

            # Wire outer adjacency
            for i in range(n_new):
                s = slots[i]
                if poly[i][2] < 0:
                    continue
                for j in range(3):
                    a, b = tris[s].v[j], tris[s].v[(j + 1) % 3]
                    if (a == poly[i][0] and b == poly[i][1]) or \
                       (a == poly[i][1] and b == poly[i][0]):
                        tris[s].neighbors[j] = poly[i][2]
                        if poly[i][3] >= 0:
                            tris[poly[i][2]].neighbors[poly[i][3]] = s
                        break

            # Wire internal adjacency
            pidx_edges = {}
            for i in range(n_new):
                s = slots[i]
                for j in range(3):
                    a, b = tris[s].v[j], tris[s].v[(j + 1) % 3]
                    if a != pidx and b != pidx:
                        continue
                    other = b if a == pidx else a
                    if other in pidx_edges:
                        os, oe = pidx_edges[other]
                        tris[s].neighbors[j] = os
                        tris[os].neighbors[oe] = s
                        del pidx_edges[other]
                    else:
                        pidx_edges[other] = (s, j)


# ============================================================
# CDT: enforce PSLG segments via cavity retriangulation
# ============================================================

def enforce_constraints(mesh, quiet):
    if not mesh.segments:
        return 0

    def build_v2t():
        v2t = [[] for _ in range(len(mesh.vertices))]
        for i, t in enumerate(mesh.triangles):
            if t.v[0] < 0:
                continue
            for j in range(3):
                v2t[t.v[j]].append(i)
        return v2t

    def ear_clip(poly_verts):
        result = []
        p = list(poly_verts)
        while len(p) > 3:
            clipped = False
            nn = len(p)
            for pass2 in range(2):
                if clipped:
                    break
                tol = 0.0 if pass2 == 0 else -1e-8
                for i in range(nn):
                    prev_v = p[(i + nn - 1) % nn]
                    cur_v = p[i]
                    next_v = p[(i + 1) % nn]
                    if orient2d(mesh.vertices[prev_v], mesh.vertices[cur_v],
                                mesh.vertices[next_v]) <= tol:
                        continue
                    inside = False
                    for k in range(nn):
                        if k in ((i + nn - 1) % nn, i, (i + 1) % nn):
                            continue
                        vi = p[k]
                        if (orient2d(mesh.vertices[prev_v], mesh.vertices[cur_v],
                                     mesh.vertices[vi]) > 0 and
                            orient2d(mesh.vertices[cur_v], mesh.vertices[next_v],
                                     mesh.vertices[vi]) > 0 and
                            orient2d(mesh.vertices[next_v], mesh.vertices[prev_v],
                                     mesh.vertices[vi]) > 0):
                            inside = True
                            break
                    if inside:
                        continue
                    result.append([prev_v, cur_v, next_v])
                    p.pop(i)
                    clipped = True
                    break
            if not clipped:
                break
        if len(p) >= 3:
            if len(p) == 3:
                o = orient2d(mesh.vertices[p[0]], mesh.vertices[p[1]], mesh.vertices[p[2]])
                if o > 0:
                    result.append([p[0], p[1], p[2]])
                elif o < 0:
                    result.append([p[0], p[2], p[1]])
            else:
                for i in range(1, len(p) - 1):
                    o = orient2d(mesh.vertices[p[0]], mesh.vertices[p[i]], mesh.vertices[p[i + 1]])
                    if o > 0:
                        result.append([p[0], p[i], p[i + 1]])
                    elif o < 0:
                        result.append([p[0], p[i + 1], p[i]])
        return result

    v2t = build_v2t()

    # Split segments at collinear intermediate vertices
    # Sort vertex indices by x-coordinate for efficient range queries
    xorder = sorted(range(len(mesh.vertices)), key=lambda i: mesh.vertices[i].x)
    xsorted = [mesh.vertices[i].x for i in xorder]

    new_segs = []
    for seg in mesh.segments:
        su, sv = seg.v0, seg.v1
        sx, sy = mesh.vertices[su].x, mesh.vertices[su].y
        ddx, ddy = mesh.vertices[sv].x - sx, mesh.vertices[sv].y - sy
        len2 = ddx * ddx + ddy * ddy
        if len2 < EPS * EPS:
            new_segs.append(seg)
            continue
        eps2_len2 = EPS * EPS * len2  # squared collinearity threshold
        bxlo = min(sx, sx + ddx) - EPS
        bxhi = max(sx, sx + ddx) + EPS
        bylo = min(sy, sy + ddy) - EPS
        byhi = max(sy, sy + ddy) + EPS
        # Binary search for x-range [bxlo, bxhi]
        ilo = bisect.bisect_left(xsorted, bxlo)
        ihi = bisect.bisect_right(xsorted, bxhi)
        on_seg = []
        for ii in range(ilo, ihi):
            i = xorder[ii]
            if i == su or i == sv:
                continue
            vy = mesh.vertices[i].y
            if vy < bylo or vy > byhi:
                continue
            px, py = mesh.vertices[i].x - sx, vy - sy
            cross = ddx * py - ddy * px
            if cross * cross > eps2_len2:
                continue
            t_param = (ddx * px + ddy * py) / len2
            if t_param <= EPS_SCALE or t_param >= 1.0 - EPS_SCALE:
                continue
            on_seg.append((t_param, i))
        if not on_seg:
            new_segs.append(seg)
        else:
            on_seg.sort()
            prev_v = su
            for _, vi in on_seg:
                new_segs.append(Segment(prev_v, vi, seg.marker))
                prev_v = vi
            new_segs.append(Segment(prev_v, sv, seg.marker))
    if len(new_segs) != len(mesh.segments):
        if not quiet:
            print(f"  Split {len(mesh.segments)} segments into {len(new_segs)} sub-segments (collinear vertices)",
                  file=sys.stderr)
        mesh.segments = new_segs

    for pass_num in range(5):
        enforced = 0
        for seg in mesh.segments:
            u, v = seg.v0, seg.v1
            if mesh.has_edge(u, v, v2t):
                continue

            success = False
            for attempt in range(4):
                if success:
                    break
                su = u if attempt % 2 == 0 else v
                sv = v if attempt % 2 == 0 else u

                sv_pert = Point(mesh.vertices[sv].x, mesh.vertices[sv].y)
                if attempt >= 2:
                    ddx = mesh.vertices[sv].x - mesh.vertices[su].x
                    ddy = mesh.vertices[sv].y - mesh.vertices[su].y
                    seg_len = math.sqrt(ddx * ddx + ddy * ddy)
                    if seg_len > EPS:
                        eps_pert = seg_len * 1e-8
                        sv_pert.y += eps_pert

                def orient_seg(p):
                    return orient2d(mesh.vertices[su], sv_pert, p)

                # Find starting triangle
                start_tri = -1
                cur_tri = -1
                for ti in v2t[su]:
                    if mesh.triangles[ti].v[0] >= 0:
                        cur_tri = ti
                        break
                if cur_tri < 0:
                    continue

                # Check if sv is already adjacent
                for ti in v2t[su]:
                    t = mesh.triangles[ti]
                    if t.v[0] < 0:
                        continue
                    for j in range(3):
                        if t.v[j] == sv:
                            success = True
                            break
                    if success:
                        break
                if success:
                    continue

                # Scan fan triangles
                best_tri = -1
                best_dot = -1e30
                sdx = mesh.vertices[sv].x - mesh.vertices[su].x
                sdy = mesh.vertices[sv].y - mesh.vertices[su].y
                for ti in v2t[su]:
                    t = mesh.triangles[ti]
                    if t.v[0] < 0:
                        continue
                    lu = -1
                    for j in range(3):
                        if t.v[j] == su:
                            lu = j
                            break
                    if lu < 0:
                        continue
                    right_v = t.v[(lu + 1) % 3]
                    left_v = t.v[(lu + 2) % 3]
                    o_right = orient_seg(mesh.vertices[right_v])
                    o_left = orient_seg(mesh.vertices[left_v])
                    if o_right >= 0 and o_left <= 0:
                        mx = (mesh.vertices[right_v].x + mesh.vertices[left_v].x) / 2
                        my = (mesh.vertices[right_v].y + mesh.vertices[left_v].y) / 2
                        dmx = mx - mesh.vertices[su].x
                        dmy = my - mesh.vertices[su].y
                        dot = sdx * dmx + sdy * dmy
                        if dot > best_dot:
                            best_dot = dot
                            best_tri = ti
                start_tri = best_tri
                if success:
                    continue
                if start_tri < 0:
                    continue

                # Walk from su to sv
                st = mesh.triangles[start_tri]
                local_u = -1
                for j in range(3):
                    if st.v[j] == su:
                        local_u = j
                a_v = st.v[(local_u + 1) % 3]
                b_v = st.v[(local_u + 2) % 3]
                oa = orient_seg(mesh.vertices[a_v])

                if oa >= 0:
                    first_above, first_below = a_v, b_v
                else:
                    first_above, first_below = b_v, a_v

                crossed_tris = [start_tri]
                above_chain = [first_above]
                below_chain = [first_below]
                bnd_edges = []

                for j in range(3):
                    ev0, ev1 = st.v[j], st.v[(j + 1) % 3]
                    if (ev0 == first_above and ev1 == first_below) or \
                       (ev0 == first_below and ev1 == first_above):
                        continue
                    nb = st.neighbors[j]
                    oe = -1
                    if nb >= 0:
                        for k in range(3):
                            if mesh.triangles[nb].neighbors[k] == start_tri:
                                oe = k
                                break
                    bnd_edges.append((ev0, ev1, nb, oe))

                reached_sv = False
                for j in range(3):
                    if st.v[j] == sv:
                        reached_sv = True
                        break

                prev_above, prev_below = first_above, first_below
                visited_tris = {start_tri}

                while not reached_sv:
                    last_tri2 = crossed_tris[-1]
                    lt = mesh.triangles[last_tri2]
                    next_tri = -1
                    for j in range(3):
                        ev0, ev1 = lt.v[j], lt.v[(j + 1) % 3]
                        if (ev0 == prev_above and ev1 == prev_below) or \
                           (ev0 == prev_below and ev1 == prev_above):
                            next_tri = lt.neighbors[j]
                            break
                    if next_tri < 0 or next_tri in visited_tris:
                        break
                    visited_tris.add(next_tri)
                    crossed_tris.append(next_tri)

                    nt = mesh.triangles[next_tri]
                    new_vert = -1
                    for j in range(3):
                        if nt.v[j] != prev_above and nt.v[j] != prev_below:
                            new_vert = nt.v[j]
                    if new_vert < 0:
                        break

                    if new_vert == sv:
                        reached_sv = True
                        for j in range(3):
                            ev0, ev1 = nt.v[j], nt.v[(j + 1) % 3]
                            if (ev0 == prev_above and ev1 == prev_below) or \
                               (ev0 == prev_below and ev1 == prev_above):
                                continue
                            nb = nt.neighbors[j]
                            oe = -1
                            if nb >= 0:
                                for k in range(3):
                                    if mesh.triangles[nb].neighbors[k] == next_tri:
                                        oe = k
                                        break
                            bnd_edges.append((ev0, ev1, nb, oe))
                    else:
                        side = orient_seg(mesh.vertices[new_vert])
                        if side > 0:
                            above_chain.append(new_vert)
                            for j in range(3):
                                ev0, ev1 = nt.v[j], nt.v[(j + 1) % 3]
                                if (ev0 == prev_above and ev1 == prev_below) or \
                                   (ev0 == prev_below and ev1 == prev_above):
                                    continue
                                if (ev0 == new_vert and ev1 == prev_below) or \
                                   (ev0 == prev_below and ev1 == new_vert):
                                    continue
                                nb = nt.neighbors[j]
                                oe = -1
                                if nb >= 0:
                                    for k in range(3):
                                        if mesh.triangles[nb].neighbors[k] == next_tri:
                                            oe = k
                                            break
                                bnd_edges.append((ev0, ev1, nb, oe))
                            prev_above = new_vert
                        else:
                            below_chain.append(new_vert)
                            for j in range(3):
                                ev0, ev1 = nt.v[j], nt.v[(j + 1) % 3]
                                if (ev0 == prev_above and ev1 == prev_below) or \
                                   (ev0 == prev_below and ev1 == prev_above):
                                    continue
                                if (ev0 == prev_above and ev1 == new_vert) or \
                                   (ev0 == new_vert and ev1 == prev_above):
                                    continue
                                nb = nt.neighbors[j]
                                oe = -1
                                if nb >= 0:
                                    for k in range(3):
                                        if mesh.triangles[nb].neighbors[k] == next_tri:
                                            oe = k
                                            break
                                bnd_edges.append((ev0, ev1, nb, oe))
                            prev_below = new_vert

                if not reached_sv:
                    continue

                above_poly = [su] + above_chain + [sv]
                below_poly = [sv] + below_chain[::-1] + [su]

                above_tris = ear_clip(above_poly)
                below_tris = ear_clip(below_poly)

                total_new = len(above_tris) + len(below_tris)
                total_old = len(crossed_tris)
                if total_new == 0:
                    continue

                crossed_tris.sort()
                slots = []
                for i in range(min(total_old, total_new)):
                    slots.append(crossed_tris[i])
                while len(slots) < total_new:
                    slots.append(len(mesh.triangles))
                    mesh.triangles.append(Triangle())
                for i in range(total_new, total_old):
                    mesh.triangles[crossed_tris[i]].v = [-1, -1, -1]
                    mesh.triangles[crossed_tris[i]].neighbors = [-1, -1, -1]

                all_new = above_tris + below_tris
                for i in range(total_new):
                    s = slots[i]
                    mesh.triangles[s].v = list(all_new[i])
                    mesh.triangles[s].neighbors = [-1, -1, -1]
                    mesh.triangles[s].region_attrib = 0
                    mesh.triangles[s].region_max_area = -1
                    mesh.triangles[s].marker = 0

                # Wire internal adjacency
                edge_map = {}
                for i in range(total_new):
                    s = slots[i]
                    t = mesh.triangles[s]
                    for j in range(3):
                        key = edge_key(t.v[j], t.v[(j + 1) % 3])
                        if key in edge_map:
                            os, oe = edge_map[key]
                            t.neighbors[j] = os
                            mesh.triangles[os].neighbors[oe] = s
                            del edge_map[key]
                        else:
                            edge_map[key] = (s, j)

                # Wire boundary edges
                for be in bnd_edges:
                    key = edge_key(be[0], be[1])
                    if key in edge_map:
                        s, j = edge_map[key]
                        mesh.triangles[s].neighbors[j] = be[2]
                        if be[2] >= 0 and be[3] >= 0:
                            mesh.triangles[be[2]].neighbors[be[3]] = s

                v2t = build_v2t()
                success = True
                enforced += 1

        if enforced == 0 and pass_num > 0:
            break

    # Direct adjacency check: flip-based
    mesh.rebuild_adjacency()
    v2t = build_v2t()
    any_fixed = True
    while any_fixed:
        any_fixed = False
        for seg in mesh.segments:
            su2, sv2 = seg.v0, seg.v1
            if mesh.has_edge(su2, sv2, v2t):
                continue
            found = False
            for ti in v2t[su2]:
                if mesh.has_edge(su2, sv2, v2t):
                    break
                t = mesh.triangles[ti]
                if t.v[0] < 0:
                    continue
                for j in range(3):
                    nb = t.neighbors[j]
                    if nb < 0:
                        continue
                    tn = mesh.triangles[nb]
                    has_sv = any(tn.v[k] == sv2 for k in range(3))
                    if not has_sv:
                        continue
                    a, b = t.v[j], t.v[(j + 1) % 3]
                    if a == su2 or b == su2 or a == sv2 or b == sv2:
                        continue
                    p_v, q_v = su2, sv2
                    o1 = orient2d(mesh.vertices[p_v], mesh.vertices[a], mesh.vertices[q_v])
                    o2 = orient2d(mesh.vertices[p_v], mesh.vertices[q_v], mesh.vertices[b])
                    if o1 < 0 or o2 < 0:
                        continue
                    # Do the flip
                    n_pa, n_bp = -1, -1
                    for k in range(3):
                        e0, e1 = t.v[k], t.v[(k + 1) % 3]
                        if (e0 == p_v and e1 == a) or (e0 == a and e1 == p_v):
                            n_pa = t.neighbors[k]
                        if (e0 == b and e1 == p_v) or (e0 == p_v and e1 == b):
                            n_bp = t.neighbors[k]
                    n_aq, n_qb = -1, -1
                    for k in range(3):
                        e0, e1 = tn.v[k], tn.v[(k + 1) % 3]
                        if (e0 == a and e1 == q_v) or (e0 == q_v and e1 == a):
                            n_aq = tn.neighbors[k]
                        if (e0 == q_v and e1 == b) or (e0 == b and e1 == q_v):
                            n_qb = tn.neighbors[k]
                    t.v = [p_v, a, q_v]
                    t.neighbors = [n_pa, n_aq, nb]
                    tn.v = [p_v, q_v, b]
                    tn.neighbors = [ti, n_qb, n_bp]
                    if n_aq >= 0:
                        for k in range(3):
                            if mesh.triangles[n_aq].neighbors[k] == nb:
                                mesh.triangles[n_aq].neighbors[k] = ti
                                break
                    if n_bp >= 0:
                        for k in range(3):
                            if mesh.triangles[n_bp].neighbors[k] == ti:
                                mesh.triangles[n_bp].neighbors[k] = nb
                                break
                    v2t = build_v2t()
                    any_fixed = True
                    found = True
                    break
                if found:
                    break

    # Edge flipping and cavity fallback for remaining segments
    seg_edges = set()
    for s2 in mesh.segments:
        if mesh.has_edge(s2.v0, s2.v1, v2t):
            seg_edges.add(edge_key(s2.v0, s2.v1))

    def segs_cross(v1, v2, v3, v4):
        d1 = orient2d(mesh.vertices[v1], mesh.vertices[v2], mesh.vertices[v3])
        d2 = orient2d(mesh.vertices[v1], mesh.vertices[v2], mesh.vertices[v4])
        if d1 * d2 >= 0:
            return False
        d3 = orient2d(mesh.vertices[v3], mesh.vertices[v4], mesh.vertices[v1])
        d4 = orient2d(mesh.vertices[v3], mesh.vertices[v4], mesh.vertices[v2])
        if d3 * d4 >= 0:
            return False
        return True

    for seg in mesh.segments:
        su, sv = seg.v0, seg.v1
        if mesh.has_edge(su, sv, v2t):
            continue

        # Walk from su toward sv collecting crossing edges, then flip them
        for attempt in range(2):
            if mesh.has_edge(su, sv, v2t):
                break
            wu = su if attempt == 0 else sv
            wv = sv if attempt == 0 else su

            # Find starting triangle in wu's fan that straddles the segment
            sdx = mesh.vertices[wv].x - mesh.vertices[wu].x
            sdy = mesh.vertices[wv].y - mesh.vertices[wu].y
            start_tri = -1
            for ti in v2t[wu]:
                t = mesh.triangles[ti]
                if t.v[0] < 0:
                    continue
                lu = -1
                for j in range(3):
                    if t.v[j] == wu:
                        lu = j
                        break
                if lu is None or lu < 0:
                    continue
                right, left = t.v[(lu + 1) % 3], t.v[(lu + 2) % 3]
                o_r = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[right])
                o_l = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[left])
                if o_r >= 0 and o_l <= 0:
                    mx = (mesh.vertices[right].x + mesh.vertices[left].x) * 0.5
                    my = (mesh.vertices[right].y + mesh.vertices[left].y) * 0.5
                    if sdx * (mx - mesh.vertices[wu].x) + sdy * (my - mesh.vertices[wu].y) > 0:
                        start_tri = ti
                        break
            if start_tri < 0:
                continue

            # Walk to collect crossing edges
            st = mesh.triangles[start_tri]
            lu = -1
            for j in range(3):
                if st.v[j] == wu:
                    lu = j
                    break
            ra, rb = st.v[(lu + 1) % 3], st.v[(lu + 2) % 3]
            oa = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[ra])
            prev_above = ra if oa >= 0 else rb
            prev_below = rb if oa >= 0 else ra

            flip_queue = deque()
            flip_queue.append(edge_key(prev_above, prev_below))

            visited = {start_tri}
            cur_tri = start_tri
            reached = False
            for _step in range(len(mesh.triangles) * 2):
                if reached:
                    break
                next_tri = -1
                ct = mesh.triangles[cur_tri]
                for j in range(3):
                    ev0, ev1 = ct.v[j], ct.v[(j + 1) % 3]
                    if (ev0 == prev_above and ev1 == prev_below) or \
                       (ev0 == prev_below and ev1 == prev_above):
                        next_tri = ct.neighbors[j]
                        break
                if next_tri < 0 or next_tri in visited:
                    break
                visited.add(next_tri)
                cur_tri = next_tri

                nt = mesh.triangles[next_tri]
                new_vert = -1
                for j in range(3):
                    if nt.v[j] != prev_above and nt.v[j] != prev_below:
                        new_vert = nt.v[j]
                        break
                if new_vert < 0:
                    break
                if new_vert == wv:
                    reached = True
                    break

                side = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[new_vert])
                if side > 0:
                    flip_queue.append(edge_key(new_vert, prev_below))
                    prev_above = new_vert
                else:
                    flip_queue.append(edge_key(prev_above, new_vert))
                    prev_below = new_vert
            if not reached:
                continue

            # Flip crossing edges (standard CDT flip algorithm)
            max_iter = len(flip_queue) * len(flip_queue) + 10
            itr = 0
            while flip_queue and itr < max_iter:
                itr += 1
                ea, eb = flip_queue.popleft()
                if not segs_cross(su, sv, ea, eb):
                    continue
                if edge_key(ea, eb) in seg_edges:
                    continue

                # Find two triangles sharing edge (ea, eb) via v2t
                ti_found, nb_found = -1, -1
                ti_edge, nb_edge2 = -1, -1
                for t_idx in v2t[ea]:
                    tr = mesh.triangles[t_idx]
                    if tr.v[0] < 0:
                        continue
                    for j in range(3):
                        if (tr.v[j] == ea and tr.v[(j + 1) % 3] == eb) or \
                           (tr.v[j] == eb and tr.v[(j + 1) % 3] == ea):
                            if ti_found < 0:
                                ti_found, ti_edge = t_idx, j
                            else:
                                nb_found, nb_edge2 = t_idx, j
                if ti_found < 0 or nb_found < 0:
                    continue

                p_v = mesh.triangles[ti_found].v[(ti_edge + 2) % 3]
                q_v = mesh.triangles[nb_found].v[(nb_edge2 + 2) % 3]

                # Check convexity
                co1 = orient2d(mesh.vertices[p_v], mesh.vertices[ea], mesh.vertices[q_v])
                co2 = orient2d(mesh.vertices[p_v], mesh.vertices[q_v], mesh.vertices[eb])
                if co1 < 0 or co2 < 0:
                    flip_queue.append((ea, eb))  # defer
                    continue

                # Flip (ea, eb) -> (p_v, q_v)
                T1, T2 = mesh.triangles[ti_found], mesh.triangles[nb_found]
                n_pa, n_bp = -1, -1
                for k in range(3):
                    e0, e1 = T1.v[k], T1.v[(k + 1) % 3]
                    if (e0 == p_v and e1 == ea) or (e0 == ea and e1 == p_v):
                        n_pa = T1.neighbors[k]
                    if (e0 == eb and e1 == p_v) or (e0 == p_v and e1 == eb):
                        n_bp = T1.neighbors[k]
                n_aq, n_qb = -1, -1
                for k in range(3):
                    e0, e1 = T2.v[k], T2.v[(k + 1) % 3]
                    if (e0 == ea and e1 == q_v) or (e0 == q_v and e1 == ea):
                        n_aq = T2.neighbors[k]
                    if (e0 == q_v and e1 == eb) or (e0 == eb and e1 == q_v):
                        n_qb = T2.neighbors[k]
                T1.v = [p_v, ea, q_v]
                T1.neighbors = [n_pa, n_aq, nb_found]
                T2.v = [p_v, q_v, eb]
                T2.neighbors = [ti_found, n_qb, n_bp]
                if n_aq >= 0:
                    for k in range(3):
                        if mesh.triangles[n_aq].neighbors[k] == nb_found:
                            mesh.triangles[n_aq].neighbors[k] = ti_found
                            break
                if n_bp >= 0:
                    for k in range(3):
                        if mesh.triangles[n_bp].neighbors[k] == ti_found:
                            mesh.triangles[n_bp].neighbors[k] = nb_found
                            break

                # Incremental v2t: T1 lost eb gained q_v, T2 lost ea gained p_v
                if ti_found in v2t[eb]:
                    v2t[eb].remove(ti_found)
                if nb_found in v2t[ea]:
                    v2t[ea].remove(nb_found)
                v2t[q_v].append(ti_found)
                v2t[p_v].append(nb_found)

                # If new edge still crosses, re-enqueue
                if segs_cross(su, sv, p_v, q_v):
                    flip_queue.append(edge_key(p_v, q_v))

        if mesh.has_edge(su, sv, v2t):
            seg_edges.add(edge_key(su, sv))
            continue

        # Cavity-based fallback
        crossed_set = set()
        for ti in range(len(mesh.triangles)):
            t = mesh.triangles[ti]
            if t.v[0] < 0:
                continue
            for j in range(3):
                a, b = t.v[j], t.v[(j + 1) % 3]
                if segs_cross(su, sv, a, b):
                    crossed_set.add(ti)
                    break
        for vi in (su, sv):
            for ti in v2t[vi]:
                t = mesh.triangles[ti]
                if t.v[0] < 0:
                    continue
                for j in range(3):
                    if t.neighbors[j] in crossed_set:
                        for k in range(3):
                            a, b = t.v[k], t.v[(k + 1) % 3]
                            if segs_cross(su, sv, a, b):
                                crossed_set.add(ti)
                                break
                        break

        if crossed_set:
            bnd_edges_cav = []
            for ti in crossed_set:
                t = mesh.triangles[ti]
                for j in range(3):
                    nb = t.neighbors[j]
                    if nb < 0 or nb not in crossed_set:
                        bnd_edges_cav.append((t.v[j], t.v[(j + 1) % 3], nb))

            # Chain boundary edges
            adj_map = {}
            for i, be in enumerate(bnd_edges_cav):
                adj_map.setdefault(be[0], []).append((be[1], i))
            for v_key in adj_map:
                if len(adj_map[v_key]) > 1:
                    adj_map[v_key].sort(
                        key=lambda x: math.atan2(
                            mesh.vertices[x[0]].y - mesh.vertices[v_key].y,
                            mesh.vertices[x[0]].x - mesh.vertices[v_key].x))

            poly_chain = []
            used_edges = set()
            start_v = bnd_edges_cav[0][0]
            cur_v = start_v
            prev_v = -1
            chain_ok = True

            for _ in range(len(bnd_edges_cav) + 1):
                poly_chain.append(cur_v)
                outs = adj_map.get(cur_v, [])
                best_idx = -1
                if len(outs) == 1 and outs[0][1] not in used_edges:
                    best_idx = 0
                elif prev_v >= 0:
                    in_dx = mesh.vertices[cur_v].x - mesh.vertices[prev_v].x
                    in_dy = mesh.vertices[cur_v].y - mesh.vertices[prev_v].y
                    in_angle = math.atan2(in_dy, in_dx)
                    best_turn = 1e30
                    for k in range(len(outs)):
                        if outs[k][1] in used_edges:
                            continue
                        ox = mesh.vertices[outs[k][0]].x - mesh.vertices[cur_v].x
                        oy = mesh.vertices[outs[k][0]].y - mesh.vertices[cur_v].y
                        out_angle = math.atan2(oy, ox)
                        turn = out_angle - in_angle
                        while turn <= 0:
                            turn += 2 * math.pi
                        if turn < best_turn:
                            best_turn = turn
                            best_idx = k
                else:
                    for k in range(len(outs)):
                        if outs[k][1] not in used_edges:
                            best_idx = k
                            break
                if best_idx < 0:
                    chain_ok = False
                    break
                used_edges.add(outs[best_idx][1])
                prev_v = cur_v
                cur_v = outs[best_idx][0]
                if cur_v == start_v and len(poly_chain) > 2:
                    break

            if len(poly_chain) >= 3 and cur_v != start_v:
                chain_ok = False

            new_tris_cav = []
            if chain_ok and len(poly_chain) >= 3:
                su_idx = sv_idx = -1
                for i in range(len(poly_chain)):
                    if poly_chain[i] == su:
                        su_idx = i
                    if poly_chain[i] == sv:
                        sv_idx = i
                if su_idx >= 0 and sv_idx >= 0:
                    nn = len(poly_chain)
                    poly1, poly2 = [], []
                    i = su_idx
                    while True:
                        poly1.append(poly_chain[i])
                        if i == sv_idx:
                            break
                        i = (i + 1) % nn
                    i = sv_idx
                    while True:
                        poly2.append(poly_chain[i])
                        if i == su_idx:
                            break
                        i = (i + 1) % nn
                    new_tris_cav = ear_clip(poly1) + ear_clip(poly2)
                else:
                    new_tris_cav = ear_clip(poly_chain)

            if new_tris_cav:
                crossed_vec = sorted(crossed_set)
                total_new = len(new_tris_cav)
                total_old = len(crossed_vec)

                slots = []
                for i in range(min(total_old, total_new)):
                    slots.append(crossed_vec[i])
                while len(slots) < total_new:
                    slots.append(len(mesh.triangles))
                    mesh.triangles.append(Triangle())
                for i in range(total_new, total_old):
                    mesh.triangles[crossed_vec[i]].v = [-1, -1, -1]
                    mesh.triangles[crossed_vec[i]].neighbors = [-1, -1, -1]

                for i in range(total_new):
                    s = slots[i]
                    mesh.triangles[s].v = list(new_tris_cav[i])
                    mesh.triangles[s].neighbors = [-1, -1, -1]

                edge_map2 = {}
                for i in range(total_new):
                    s = slots[i]
                    t = mesh.triangles[s]
                    for j in range(3):
                        key = edge_key(t.v[j], t.v[(j + 1) % 3])
                        if key in edge_map2:
                            os, oe = edge_map2[key]
                            t.neighbors[j] = os
                            mesh.triangles[os].neighbors[oe] = s
                            del edge_map2[key]
                        else:
                            edge_map2[key] = (s, j)

                for be in bnd_edges_cav:
                    key = edge_key(be[0], be[1])
                    if key in edge_map2:
                        s, j = edge_map2[key]
                        mesh.triangles[s].neighbors[j] = be[2]
                        if be[2] >= 0:
                            ot = mesh.triangles[be[2]]
                            for k in range(3):
                                if ot.neighbors[k] in crossed_set:
                                    e0, e1 = ot.v[k], ot.v[(k + 1) % 3]
                                    if edge_key(e0, e1) == key:
                                        ot.neighbors[k] = s
                                        break
                v2t = build_v2t()

        if mesh.has_edge(su, sv, v2t):
            seg_edges.add(edge_key(su, sv))

    # Count and return remaining missing segments
    miss = 0
    v2t = build_v2t()
    for s in mesh.segments:
        if not mesh.has_edge(s.v0, s.v1, v2t):
            miss += 1
            if not quiet:
                print(f"  UNENFORCED: v{s.v0} ({mesh.vertices[s.v0].x},{mesh.vertices[s.v0].y})"
                      f"->v{s.v1} ({mesh.vertices[s.v1].x},{mesh.vertices[s.v1].y})"
                      f" len={math.sqrt((mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)*(mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)+(mesh.vertices[s.v1].y-mesh.vertices[s.v0].y)*(mesh.vertices[s.v1].y-mesh.vertices[s.v0].y))}",
                      file=sys.stderr)
    if not quiet and miss > 0:
        print(f"Warning: {miss} of {len(mesh.segments)} segments could not be enforced.",
              file=sys.stderr)
    return miss


# ============================================================
# Flood-fill hole removal
# ============================================================

def remove_holes(mesh, opts):
    if opts.no_holes:
        return
    mesh.rebuild_adjacency()

    seg_edges = set()
    for s in mesh.segments:
        seg_edges.add(edge_key(s.v0, s.v1))

    remove = [False] * len(mesh.triangles)

    # Remove exterior: flood-fill from boundary triangles whose boundary
    # edge is not a constrained segment.
    q = deque()
    for i in range(len(mesh.triangles)):
        if remove[i]:
            continue
        t = mesh.triangles[i]
        for j in range(3):
            if t.neighbors[j] >= 0:
                continue
            key = edge_key(t.v[j], t.v[(j + 1) % 3])
            if key in seg_edges:
                continue
            if not remove[i]:
                remove[i] = True
                q.append(i)
            break
    while q:
        ti = q.popleft()
        t = mesh.triangles[ti]
        for j in range(3):
            nb = t.neighbors[j]
            if nb < 0 or remove[nb]:
                continue
            key = edge_key(t.v[j], t.v[(j + 1) % 3])
            if key in seg_edges:
                continue
            remove[nb] = True
            q.append(nb)

    # Remove explicit holes via flood fill from seed points
    for hole in mesh.holes:
        seed_tri = mesh.locate_triangle(hole.x, hole.y)
        if seed_tri < 0:
            continue
        q = deque([seed_tri])
        remove[seed_tri] = True
        while q:
            ti = q.popleft()
            t = mesh.triangles[ti]
            for j in range(3):
                nb = t.neighbors[j]
                if nb < 0 or remove[nb]:
                    continue
                key = edge_key(t.v[j], t.v[(j + 1) % 3])
                if key in seg_edges:
                    continue
                remove[nb] = True
                q.append(nb)

    kept = []
    remap = [-1] * len(mesh.triangles)
    for i in range(len(mesh.triangles)):
        if not remove[i]:
            remap[i] = len(kept)
            kept.append(mesh.triangles[i])
    for t in kept:
        for j in range(3):
            t.neighbors[j] = remap[t.neighbors[j]] if t.neighbors[j] >= 0 else -1
    mesh.triangles = kept


# ============================================================
# Flood-fill region assignment
# ============================================================

def assign_regions(mesh, opts):
    if not opts.regions or not mesh.regions:
        return
    mesh.rebuild_adjacency()

    seg_edges = set()
    for s in mesh.segments:
        seg_edges.add(edge_key(s.v0, s.v1))

    for t in mesh.triangles:
        t.region_attrib = 0
        t.region_max_area = -1

    for reg in mesh.regions:
        seed_tri = mesh.locate_triangle(reg.x, reg.y)
        if seed_tri < 0:
            continue
        q = deque([seed_tri])
        visited = [False] * len(mesh.triangles)
        visited[seed_tri] = True
        while q:
            ti = q.popleft()
            t = mesh.triangles[ti]
            t.region_attrib = reg.attrib
            t.region_max_area = reg.max_area
            for j in range(3):
                nb = t.neighbors[j]
                if nb < 0 or visited[nb]:
                    continue
                key = edge_key(t.v[j], t.v[(j + 1) % 3])
                if key in seg_edges:
                    continue
                visited[nb] = True
                q.append(nb)


# ============================================================
# Point insertion for quality refinement (Lawson flip-based)
# ============================================================

def insert_point_lawson(mesh, np, hint_tri, constrained_edges, affected):
    pidx = len(mesh.vertices)
    mesh.vertices.append(np)

    ti = mesh.locate_triangle(np.x, np.y, hint_tri)
    if ti < 0:
        mesh.vertices.pop()
        return -1

    t = mesh.triangles[ti]
    for j in range(3):
        dx = mesh.vertices[t.v[j]].x - np.x
        dy = mesh.vertices[t.v[j]].y - np.y
        if dx * dx + dy * dy < EPS * EPS:
            mesh.vertices.pop()
            return -1

    # Inherit region from the containing triangle
    reg_attr = t.region_attrib
    reg_area = t.region_max_area

    # Check if point is on an edge
    on_edge = -1
    for j in range(3):
        va, vb = t.v[j], t.v[(j + 1) % 3]
        edx = mesh.vertices[vb].x - mesh.vertices[va].x
        edy = mesh.vertices[vb].y - mesh.vertices[va].y
        o = orient2d_xy(mesh.vertices[va].x, mesh.vertices[va].y,
                        mesh.vertices[vb].x, mesh.vertices[vb].y, np.x, np.y)
        if o * o < EPS * EPS * (edx * edx + edy * edy):
            on_edge = j
            break

    flip_stack = []

    if on_edge < 0:
        # Interior: split triangle into 3
        v0, v1, v2 = t.v[0], t.v[1], t.v[2]
        n0, n1, n2 = t.neighbors[0], t.neighbors[1], t.neighbors[2]

        t0 = ti
        t1 = len(mesh.triangles)
        mesh.triangles.append(Triangle())
        t2 = len(mesh.triangles)
        mesh.triangles.append(Triangle())

        mesh.triangles[t0] = Triangle(v=[v0, v1, pidx], neighbors=[n0, t1, t2],
                                      region_attrib=reg_attr, region_max_area=reg_area)
        mesh.triangles[t1] = Triangle(v=[v1, v2, pidx], neighbors=[n1, t2, t0],
                                      region_attrib=reg_attr, region_max_area=reg_area)
        mesh.triangles[t2] = Triangle(v=[v2, v0, pidx], neighbors=[n2, t0, t1],
                                      region_attrib=reg_attr, region_max_area=reg_area)

        if n1 >= 0:
            for k in range(3):
                if mesh.triangles[n1].neighbors[k] == ti:
                    mesh.triangles[n1].neighbors[k] = t1
                    break
        if n2 >= 0:
            for k in range(3):
                if mesh.triangles[n2].neighbors[k] == ti:
                    mesh.triangles[n2].neighbors[k] = t2
                    break

        flip_stack.extend([(t0, 0), (t1, 0), (t2, 0)])
        affected.extend([t0, t1, t2])
    else:
        # On edge: split 2 triangles into 4
        va = t.v[on_edge]
        vb = t.v[(on_edge + 1) % 3]
        vc = t.v[(on_edge + 2) % 3]
        n_ab = t.neighbors[on_edge]
        n_bvc = t.neighbors[(on_edge + 1) % 3]
        n_cva = t.neighbors[(on_edge + 2) % 3]

        if n_ab < 0:
            # Boundary edge: split 1 triangle into 2
            t0 = ti
            t1 = len(mesh.triangles)
            mesh.triangles.append(Triangle())

            mesh.triangles[t0] = Triangle(v=[va, pidx, vc], neighbors=[-1, t1, n_cva],
                                          region_attrib=reg_attr, region_max_area=reg_area)
            mesh.triangles[t1] = Triangle(v=[pidx, vb, vc], neighbors=[-1, n_bvc, t0],
                                          region_attrib=reg_attr, region_max_area=reg_area)

            if n_bvc >= 0:
                for k in range(3):
                    if mesh.triangles[n_bvc].neighbors[k] == ti:
                        mesh.triangles[n_bvc].neighbors[k] = t1
                        break

            flip_stack.extend([(t0, 2), (t1, 1)])
            affected.extend([t0, t1])
        else:
            # Both sides: split 2 triangles into 4
            tn = mesh.triangles[n_ab]
            je = -1
            for k in range(3):
                if tn.neighbors[k] == ti:
                    je = k
                    break
            if je < 0:
                mesh.vertices.pop()
                return -1

            vd = tn.v[(je + 2) % 3]
            if tn.v[je] == vb and tn.v[(je + 1) % 3] == va:
                n_avd = tn.neighbors[(je + 1) % 3]
                n_dvb = tn.neighbors[(je + 2) % 3]
            else:
                n_dvb = tn.neighbors[(je + 1) % 3]
                n_avd = tn.neighbors[(je + 2) % 3]

            # Inherit region from the neighbor triangle
            reg_attr2 = tn.region_attrib
            reg_area2 = tn.region_max_area

            t0 = ti
            t1 = len(mesh.triangles)
            mesh.triangles.append(Triangle())
            t2 = n_ab
            t3 = len(mesh.triangles)
            mesh.triangles.append(Triangle())

            mesh.triangles[t0] = Triangle(v=[va, pidx, vc], neighbors=[t3, t1, n_cva],
                                          region_attrib=reg_attr, region_max_area=reg_area)
            mesh.triangles[t1] = Triangle(v=[pidx, vb, vc], neighbors=[t2, n_bvc, t0],
                                          region_attrib=reg_attr, region_max_area=reg_area)
            mesh.triangles[t2] = Triangle(v=[vb, pidx, vd], neighbors=[t1, t3, n_dvb],
                                          region_attrib=reg_attr2, region_max_area=reg_area2)
            mesh.triangles[t3] = Triangle(v=[pidx, va, vd], neighbors=[t0, n_avd, t2],
                                          region_attrib=reg_attr2, region_max_area=reg_area2)

            if n_bvc >= 0:
                for k in range(3):
                    if mesh.triangles[n_bvc].neighbors[k] == ti:
                        mesh.triangles[n_bvc].neighbors[k] = t1
                        break
            if n_avd >= 0:
                for k in range(3):
                    if mesh.triangles[n_avd].neighbors[k] == n_ab:
                        mesh.triangles[n_avd].neighbors[k] = t3
                        break

            flip_stack.extend([(t0, 2), (t1, 1), (t2, 2), (t3, 1)])
            affected.extend([t0, t1, t2, t3])

    # Lawson flip loop
    while flip_stack:
        ft, fj = flip_stack.pop()
        tri = mesh.triangles[ft]
        if tri.v[0] < 0:
            continue
        fn = tri.neighbors[fj]
        if fn < 0:
            continue

        fa, fb, fp = tri.v[fj], tri.v[(fj + 1) % 3], tri.v[(fj + 2) % 3]

        if edge_key(fa, fb) in constrained_edges:
            continue

        tri_n = mesh.triangles[fn]
        if tri_n.v[0] < 0:
            continue
        fj2 = -1
        for k in range(3):
            if tri_n.neighbors[k] == ft:
                fj2 = k
                break
        if fj2 < 0:
            continue

        fq = tri_n.v[(fj2 + 2) % 3]

        if in_circle(mesh.vertices[fa], mesh.vertices[fb], mesh.vertices[fp],
                     mesh.vertices[fq]) <= 0:
            continue

        if orient2d(mesh.vertices[fa], mesh.vertices[fq], mesh.vertices[fp]) <= 0:
            continue
        if orient2d(mesh.vertices[fb], mesh.vertices[fp], mesh.vertices[fq]) <= 0:
            continue

        n_fb_fp = tri.neighbors[(fj + 1) % 3]
        n_fp_fa = tri.neighbors[(fj + 2) % 3]
        n_ext1 = tri_n.neighbors[(fj2 + 1) % 3]
        n_ext2 = tri_n.neighbors[(fj2 + 2) % 3]

        if tri_n.v[fj2] != fb:
            n_ext1, n_ext2 = n_ext2, n_ext1

        ra, ra2 = tri.region_attrib, tri.region_max_area
        rb, rb2 = tri_n.region_attrib, tri_n.region_max_area

        mesh.triangles[ft] = Triangle(v=[fa, fq, fp], neighbors=[n_ext1, fn, n_fp_fa],
                                      region_attrib=ra, region_max_area=ra2)
        mesh.triangles[fn] = Triangle(v=[fb, fp, fq], neighbors=[n_fb_fp, ft, n_ext2],
                                      region_attrib=rb, region_max_area=rb2)

        if n_ext1 >= 0:
            for k in range(3):
                if mesh.triangles[n_ext1].neighbors[k] == fn:
                    mesh.triangles[n_ext1].neighbors[k] = ft
                    break
        if n_fb_fp >= 0:
            for k in range(3):
                if mesh.triangles[n_fb_fp].neighbors[k] == ft:
                    mesh.triangles[n_fb_fp].neighbors[k] = fn
                    break

        affected.extend([ft, fn])
        flip_stack.extend([(ft, 0), (fn, 2)])

    # Bump generation
    for ai in affected:
        if 0 <= ai < len(mesh.triangles):
            mesh.triangles[ai].generation += 1

    return pidx


# ============================================================
# Quality refinement
# ============================================================

def tri_metric(mesh, ti, min_ang_threshold, global_max_a, use_region_area,
               cos_min_ang=-1.0):
    t = mesh.triangles[ti]
    if t.v[0] < 0:
        return (0, False, False, 0, -1, 0, 0, 0, 0)
    a = mesh.vertices[t.v[0]]
    b = mesh.vertices[t.v[1]]
    c = mesh.vertices[t.v[2]]

    # Squared edge lengths (computed once, reused for both area and angle checks)
    la2 = (b.x-c.x)*(b.x-c.x) + (b.y-c.y)*(b.y-c.y)
    lb2 = (a.x-c.x)*(a.x-c.x) + (a.y-c.y)*(a.y-c.y)
    lc2 = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)

    # Area via cross product (no orient2d exact-arithmetic overhead)
    cross = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x)
    area = 0.5 * abs(cross)

    effective_max_area = -1
    if use_region_area and t.region_max_area > 0:
        effective_max_area = t.region_max_area
    if global_max_a > 0:
        if effective_max_area > 0:
            effective_max_area = min(effective_max_area, global_max_a)
        else:
            effective_max_area = global_max_a
    # LFS constraint: if any vertex has lfs defined, impose area <= sqrt(3)/4 * lfs^2
    for j in range(3):
        vl = mesh.vertices[t.v[j]].lfs
        if vl > 0:
            lfs_area = 0.43301270189221932 * vl * vl  # sqrt(3)/4
            if effective_max_area > 0:
                effective_max_area = min(effective_max_area, lfs_area)
            else:
                effective_max_area = lfs_area
    area_viol = effective_max_area > 0 and area > effective_max_area + EPS

    # Fast angle check using cosine comparison (no hypot/acos)
    # Minimum angle < threshold iff maximum cosine > cosMinAng.
    # cos(angle_opposite_la) = (lb²+lc²-la²) / (2·sqrt(lb²·lc²))
    # Compare squared: (lb²+lc²-la²)² > cos²Threshold · 4·lb²·lc²
    ang_viol = False
    max_cos_ratio = 0.0
    if min_ang_threshold > 0 and cos_min_ang > -0.5:
        cos2t = cos_min_ang * cos_min_ang
        for num, prod in ((lb2 + lc2 - la2, lb2 * lc2),
                          (la2 + lc2 - lb2, la2 * lc2),
                          (la2 + lb2 - lc2, la2 * lb2)):
            if num > 0:
                denom = cos2t * 4.0 * prod
                if denom > 0:
                    ratio = num * num / denom
                    if ratio > 1.0:
                        ang_viol = True
                        if ratio > max_cos_ratio:
                            max_cos_ratio = ratio
    elif min_ang_threshold > 0:
        ang = min_angle(a, b, c)
        ang_viol = ang < min_ang_threshold - EPS
        if ang_viol:
            max_cos_ratio = (min_ang_threshold - ang) / min_ang_threshold + 1.0

    if not area_viol and not ang_viol:
        return (0, False, False, area, effective_max_area, la2, lb2, lc2, cross)
    metric = 0
    if area_viol:
        metric = area / (effective_max_area if effective_max_area > 0 else 1.0)
    if ang_viol:
        metric = max(metric, max_cos_ratio)
    return (metric, ang_viol, area_viol, area, effective_max_area, la2, lb2, lc2, cross)


def refine_quality(mesh, opts, n_input_verts):
    if not opts.quality and not opts.area_limit:
        return
    min_ang_val = opts.min_angle + 0.1 if opts.quality else 0.0
    cos_min_ang = math.cos(min_ang_val * math.pi / 180.0) if min_ang_val > 0 else -1.0
    global_max_a = opts.max_area if opts.area_limit else -1.0
    use_region_area = opts.area_limit

    # Estimate final triangle count from CDT using LFS integral with
    # circle-of-influence model for refinement seed vertices.
    nv = len(mesh.vertices)
    nt = len(mesh.triangles)

    # Build vertex-to-triangle adjacency
    v2t = [[] for _ in range(nv)]
    for i in range(nt):
        t = mesh.triangles[i]
        if t.v[0] < 0:
            continue
        for j in range(3):
            v2t[t.v[j]].append(i)

    # Bounding box diagonal — used to cap degenerate circumradii
    gxmin = min(v.x for v in mesh.vertices)
    gxmax = max(v.x for v in mesh.vertices)
    gymin = min(v.y for v in mesh.vertices)
    gymax = max(v.y for v in mesh.vertices)
    bb_diag = math.sqrt((gxmax-gxmin)*(gxmax-gxmin) + (gymax-gymin)*(gymax-gymin))

    # Compute LFS at each vertex as min circumradius of incident CDT triangles
    lfs = [1e30] * nv
    for i in range(nt):
        t = mesh.triangles[i]
        if t.v[0] < 0:
            continue
        p0, p1, p2 = mesh.vertices[t.v[0]], mesh.vertices[t.v[1]], mesh.vertices[t.v[2]]
        area = tri_area(p0, p1, p2)
        e0 = math.sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y))
        e1 = math.sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y))
        e2 = math.sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y))
        R = e0 * e1 * e2 / (4.0 * area) if area > 0 else 0
        if R > 0:
            R = min(R, bb_diag)
            for j in range(3):
                lfs[t.v[j]] = min(lfs[t.v[j]], R)

    # Fold in user-assigned LFS constraints
    for vi in range(nv):
        if mesh.vertices[vi].lfs > 0:
            lfs[vi] = min(lfs[vi], mesh.vertices[vi].lfs)

    # Compute background mesh size H at each vertex from area constraints
    # and median circumradius of incident CDT triangles
    bg_size = [1e30] * nv
    for vi in range(nv):
        # Area constraint contribution (only if area limiting is active)
        amax = global_max_a
        if use_region_area:
            for ti in v2t[vi]:
                t = mesh.triangles[ti]
                if t.v[0] < 0:
                    continue
                if t.region_max_area > 0:
                    amax = min(amax, t.region_max_area) if amax > 0 else t.region_max_area
        if amax > 0:
            bg_size[vi] = math.sqrt(4.0 * amax / math.sqrt(3.0))

        # Median circumradius of incident CDT triangles
        inc_r = []
        for ti in v2t[vi]:
            t = mesh.triangles[ti]
            if t.v[0] < 0:
                continue
            p0, p1, p2 = mesh.vertices[t.v[0]], mesh.vertices[t.v[1]], mesh.vertices[t.v[2]]
            a2 = tri_area(p0, p1, p2)
            ee0 = math.sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y))
            ee1 = math.sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y))
            ee2 = math.sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y))
            R2 = ee0 * ee1 * ee2 / (4.0 * a2) if a2 > 0 else 0
            if R2 > 0:
                inc_r.append(min(R2, bb_diag))
        if inc_r:
            inc_r.sort()
            med_r = inc_r[len(inc_r) // 2]
            bg_size[vi] = min(bg_size[vi], med_r)

    # Pass 1: background mesh estimate from CDT triangles
    alpha = math.sqrt(3.0) / 4.0
    pi_over_alpha = math.pi / alpha
    cos_min_angle = math.cos(min_ang_val * math.pi / 180.0)

    est_background = 0.0
    n_cdt = 0
    n_bad_cdt = 0
    total_area = 0.0
    for i in range(nt):
        t = mesh.triangles[i]
        if t.v[0] < 0:
            continue
        n_cdt += 1
        p0, p1, p2 = mesh.vertices[t.v[0]], mesh.vertices[t.v[1]], mesh.vertices[t.v[2]]
        area = tri_area(p0, p1, p2)
        total_area += area
        H = (bg_size[t.v[0]] + bg_size[t.v[1]] + bg_size[t.v[2]]) / 3.0
        bg_est = area / (alpha * H * H)
        # Check for bad angle
        e0 = math.sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y))
        e1 = math.sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y))
        e2 = math.sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y))
        cos_a0 = (e1*e1 + e2*e2 - e0*e0) / (2*e1*e2 + 1e-30)
        cos_a1 = (e0*e0 + e2*e2 - e1*e1) / (2*e0*e2 + 1e-30)
        cos_a2 = (e0*e0 + e1*e1 - e2*e2) / (2*e0*e1 + 1e-30)
        cos_max = max(cos_a0, cos_a1, cos_a2)
        is_bad = cos_max > cos_min_angle
        if is_bad:
            n_bad_cdt += 1
        # Angle-area interaction: 1.7x when area dominates on bad triangles
        if is_bad and bg_est > 6.0:
            est_background += bg_est * 1.7
        else:
            est_background += max(bg_est, 6.0)

    # Pass 2: refinement seed contributions (circle-of-influence model)
    est_seeds = 0.0
    n_seeds = 0
    for vi in range(nv):
        h = lfs[vi]
        H = bg_size[vi]
        if h >= 1e29 or h >= H * 0.5:
            continue
        n_seeds += 1
        ratio = H / h
        log_r = math.log(ratio)
        # Correction: Ruppert's geometric grading creates quadratically more
        # triangles than the linear grading model at high H/h ratios.
        seed_contrib = pi_over_alpha * (1.0 + 2.0 * log_r) * max(1.0, log_r)
        seed_contrib -= pi_over_alpha  # subtract background already counted
        if seed_contrib > 0:
            est_seeds += seed_contrib

    # Pass 3: segment LFS contributions
    # A segment with LFS constraint h generates a grading band on both sides.
    # Theoretical integral (linear grading) telescopes to 2*L/(alpha*h), but
    # Ruppert's geometric grading is ~3x denser, giving 6*L/(alpha*h).
    # NOTE: this assumes grading on both sides (interior segment).  Boundary
    # segments only grade on one side, so this overestimates by ~2x for them.
    # Conservative is fine for a safety cap; if a tighter estimate is ever
    # needed, halve the contribution for boundary segments (those with only
    # one incident triangle on the constrained edge).
    est_seg_lfs = 0
    for s in mesh.segments:
        if s.lfs <= 0:
            continue
        Ldx = mesh.vertices[s.v1].x - mesh.vertices[s.v0].x
        Ldy = mesh.vertices[s.v1].y - mesh.vertices[s.v0].y
        L = math.sqrt(Ldx*Ldx + Ldy*Ldy)
        est_seg_lfs += 6.0 * L / (alpha * s.lfs)

    est_total_tris = est_background + est_seeds + est_seg_lfs
    max_steiner = max(10000, int(est_total_tris * 5.0))
    if not opts.quiet:
        sys.stderr.write(f"  CDT: {n_cdt} tris ({n_bad_cdt} bad), area={total_area:.6g}, seeds: {n_seeds}\n")
        sys.stderr.write(f"  Est triangles: {int(est_total_tris)}, Steiner limit: {max_steiner}\n")
    steiner_count = [0]  # use list for mutability in closures

    constrained_edges = set()
    for s in mesh.segments:
        constrained_edges.add(edge_key(s.v0, s.v1))

    # Detect acute apices and initialize shell tracking (concentric shells)
    acute_apices = set()
    vert_segs = {}
    for si in range(len(mesh.segments)):
        for v in (mesh.segments[si].v0, mesh.segments[si].v1):
            vert_segs.setdefault(v, []).append(si)
    for vi, slist in vert_segs.items():
        if vi >= n_input_verts or len(slist) < 2:
            continue
        for i in range(len(slist)):
            for j in range(i + 1, len(slist)):
                oi = mesh.segments[slist[i]].v1 if mesh.segments[slist[i]].v0 == vi else mesh.segments[slist[i]].v0
                oj = mesh.segments[slist[j]].v1 if mesh.segments[slist[j]].v0 == vi else mesh.segments[slist[j]].v0
                dix = mesh.vertices[oi].x - mesh.vertices[vi].x
                diy = mesh.vertices[oi].y - mesh.vertices[vi].y
                djx = mesh.vertices[oj].x - mesh.vertices[vi].x
                djy = mesh.vertices[oj].y - mesh.vertices[vi].y
                dot = dix * djx + diy * djy
                cross = dix * djy - diy * djx
                if dot > 0 and abs(cross) < dot * 1.15:  # angle < ~60°
                    acute_apices.add(vi)
    shell_apex = [-1] * len(mesh.vertices)

    gxmin -= 1
    gymin -= 1
    gxmax += 1
    gymax += 1
    g_nx, g_ny = 50, 50
    gcell = max((gxmax - gxmin) / g_nx, (gymax - gymin) / g_ny)
    g_nx = int((gxmax - gxmin) / gcell) + 1
    g_ny = int((gymax - gymin) / gcell) + 1

    seg_grid = {}

    def grid_add_seg(si):
        s = mesh.segments[si]
        sa, sb = mesh.vertices[s.v0], mesh.vertices[s.v1]
        mx, my = (sa.x + sb.x) * 0.5, (sa.y + sb.y) * 0.5
        r = math.sqrt((sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y)) * 0.5
        cx0 = max(0, int((mx - r - gxmin) / gcell))
        cy0 = max(0, int((my - r - gymin) / gcell))
        cx1 = min(g_nx - 1, int((mx + r - gxmin) / gcell))
        cy1 = min(g_ny - 1, int((my + r - gymin) / gcell))
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                key = cx * g_ny + cy
                seg_grid.setdefault(key, []).append(si)

    def grid_remove_seg(si):
        s = mesh.segments[si]
        sa, sb = mesh.vertices[s.v0], mesh.vertices[s.v1]
        mx, my = (sa.x + sb.x) * 0.5, (sa.y + sb.y) * 0.5
        r = math.sqrt((sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y)) * 0.5
        cx0 = max(0, int((mx - r - gxmin) / gcell))
        cy0 = max(0, int((my - r - gymin) / gcell))
        cx1 = min(g_nx - 1, int((mx + r - gxmin) / gcell))
        cy1 = min(g_ny - 1, int((my + r - gymin) / gcell))
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                key = cx * g_ny + cy
                if key in seg_grid:
                    lst = seg_grid[key]
                    if si in lst:
                        lst.remove(si)

    # edge_to_seg map: edge_key -> segment index (for PBC partner lookup)
    edge_to_seg = {}
    for si in range(len(mesh.segments)):
        grid_add_seg(si)
        edge_to_seg[edge_key(mesh.segments[si].v0, mesh.segments[si].v1)] = si

    def split_segment_core(si, out_affected, hint_tri=-1):
        """Split segment at midpoint. Returns midpoint index or -1 on failure."""
        seg = Segment(mesh.segments[si].v0, mesh.segments[si].v1,
                      mesh.segments[si].marker, mesh.segments[si].lfs,
                      mesh.segments[si].pbc_partner, mesh.segments[si].pbc_type)
        sa, sb = mesh.vertices[seg.v0], mesh.vertices[seg.v1]
        mx, my = (sa.x + sb.x) * 0.5, (sa.y + sb.y) * 0.5
        mid = Point(mx, my, len(mesh.vertices), seg.marker, [], seg.lfs)

        ek = edge_key(seg.v0, seg.v1)
        constrained_edges.discard(ek)
        edge_to_seg.pop(ek, None)
        grid_remove_seg(si)

        mid_idx = insert_point_lawson(mesh, mid, hint_tri, constrained_edges, out_affected)
        if mid_idx < 0:
            constrained_edges.add(ek)
            edge_to_seg[ek] = si
            grid_add_seg(si)
            return -1
        steiner_count[0] += 1
        # Propagate shell membership for acute-angle wedge tracking
        while len(shell_apex) <= mid_idx:
            shell_apex.append(-1)
        a0 = seg.v0 if seg.v0 < n_input_verts else shell_apex[seg.v0]
        a1 = seg.v1 if seg.v1 < n_input_verts else shell_apex[seg.v1]
        if a0 >= 0 and a0 in acute_apices:
            shell_apex[mid_idx] = a0
        elif a1 >= 0 and a1 in acute_apices:
            shell_apex[mid_idx] = a1
        s1 = Segment(seg.v0, mid_idx, seg.marker, seg.lfs, seg.pbc_partner, seg.pbc_type)
        s2 = Segment(mid_idx, seg.v1, seg.marker, seg.lfs, seg.pbc_partner, seg.pbc_type)
        mesh.segments[si] = s1
        new_si = len(mesh.segments)
        mesh.segments.append(s2)
        ek1 = edge_key(s1.v0, s1.v1)
        ek2 = edge_key(s2.v0, s2.v1)
        constrained_edges.add(ek1)
        constrained_edges.add(ek2)
        edge_to_seg[ek1] = si
        edge_to_seg[ek2] = new_si
        grid_add_seg(si)
        grid_add_seg(new_si)
        return mid_idx

    def split_segment(si, out_affected, hint_tri=-1):
        """Split segment; for PBC segments, also split the partner."""
        seg = mesh.segments[si]  # read before split modifies it
        pbc_type = seg.pbc_type
        sv0, sv1 = seg.v0, seg.v1
        mid_idx = split_segment_core(si, out_affected, hint_tri)
        if mid_idx < 0:
            return False

        # Synchronized PBC partner split
        if pbc_type >= 0:
            tv0 = mesh.pbc_twin.get(sv0, -1)
            tv1 = mesh.pbc_twin.get(sv1, -1)
            if tv0 >= 0 and tv1 >= 0:
                pek = edge_key(tv0, tv1)
                partner_si = edge_to_seg.get(pek, -1)
                if partner_si >= 0:
                    partner_affected = []
                    partner_mid = split_segment_core(partner_si, partner_affected)
                    if partner_mid >= 0:
                        mesh.pbc_twin[mid_idx] = partner_mid
                        mesh.pbc_twin[partner_mid] = mid_idx
                        mesh.pbc_node_type[mid_idx] = pbc_type
                        mesh.pbc_node_type[partner_mid] = pbc_type
        return True

    def find_encroached(px, py):
        gcx = max(0, min(g_nx - 1, int((px - gxmin) / gcell)))
        gcy = max(0, min(g_ny - 1, int((py - gymin) / gcell)))
        for ddx in range(-1, 2):
            for ddy in range(-1, 2):
                cx, cy = gcx + ddx, gcy + ddy
                if cx < 0 or cx >= g_nx or cy < 0 or cy >= g_ny:
                    continue
                key = cx * g_ny + cy
                if key not in seg_grid:
                    continue
                for si in seg_grid[key]:
                    s = mesh.segments[si]
                    sa, sb = mesh.vertices[s.v0], mesh.vertices[s.v1]
                    mx2, my2 = (sa.x + sb.x) * 0.5, (sa.y + sb.y) * 0.5
                    sr2 = ((sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y)) * 0.25
                    dd2 = (px-mx2)*(px-mx2)+(py-my2)*(py-my2)
                    if dd2 < sr2:
                        return si
        return -1

    # Priority queue: use negative metric for max-heap behavior
    pq = []
    for i in range(len(mesh.triangles)):
        mr = tri_metric(mesh, i, min_ang_val, global_max_a, use_region_area, cos_min_ang)
        if mr[0] > 0:
            heappush(pq, (-mr[0], i, mesh.triangles[i].generation))

    # Precompute off-center constant (loop-invariant)
    offconst = (0.475 * math.sqrt(
        (1.0 + math.cos(min_ang_val * math.pi / 180.0)) /
        (1.0 - math.cos(min_ang_val * math.pi / 180.0)))) if min_ang_val > 0 else 0

    dbg_encroach = dbg_insert = dbg_skip = dbg_fail = 0

    while pq and steiner_count[0] < max_steiner:
        neg_metric, bad_idx, gen = heappop(pq)

        if bad_idx >= len(mesh.triangles):
            continue
        bt = mesh.triangles[bad_idx]
        if bt.v[0] < 0:
            continue
        if bt.generation != gen:
            continue
        tmr = tri_metric(mesh, bad_idx, min_ang_val, global_max_a, use_region_area, cos_min_ang)
        if tmr[0] <= 0:
            continue
        ang_viol, area_viol, area, effective_max_area = tmr[1], tmr[2], tmr[3], tmr[4]
        tmr_la2, tmr_lb2, tmr_lc2 = tmr[5], tmr[6], tmr[7]

        a = mesh.vertices[bt.v[0]]
        b = mesh.vertices[bt.v[1]]
        c = mesh.vertices[bt.v[2]]

        # Skip if both edges at smallest-angle vertex are segments
        if ang_viol and not area_viol:
            # Find vertex with smallest angle = vertex opposite the shortest edge
            # Use cached squared edge lengths from tri_metric
            la2, lb2, lc2 = tmr_la2, tmr_lb2, tmr_lc2
            min_v = 0  # opposite la (edge b-c)
            if lb2 < la2 and lb2 < lc2:
                min_v = 1  # opposite lb (edge a-c)
            elif lc2 < la2:
                min_v = 2  # opposite lc (edge a-b)
            va = bt.v[min_v]
            vb_v = bt.v[(min_v + 1) % 3]
            vc_v = bt.v[(min_v + 2) % 3]
            if (edge_key(va, vb_v) in constrained_edges and
                edge_key(va, vc_v) in constrained_edges):
                dbg_skip += 1
                continue
            # Check if worst-angle vertex is inside a small-angle shell
            if va < len(shell_apex) and shell_apex[va] >= 0:
                dbg_skip += 1
                continue

        # Compute circumcenter or off-center
        if ang_viol:
            # Inline circumcenter using cached values from tri_metric
            xdo, ydo = b.x - a.x, b.y - a.y
            xao, yao = c.x - a.x, c.y - a.y
            D = 2.0 * (xdo * yao - ydo * xao)
            if abs(D) < EPS:
                tcx, tcy = (a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0
            else:
                # dodist=xdo²+ydo²=lc2, aodist=xao²+yao²=lb2 (from tri_metric)
                dodist_cc, aodist_cc = tmr_lc2, tmr_lb2
                tcx = a.x + (yao * dodist_cc - ydo * aodist_cc) / D
                tcy = a.y + (xdo * aodist_cc - xao * dodist_cc) / D
            dodist, aodist, dadist = tmr_lc2, tmr_lb2, tmr_la2
            ddx_cc, ddy_cc = tcx - a.x, tcy - a.y
            if dodist < aodist and dodist < dadist:
                dxoff = 0.5 * xdo - offconst * ydo
                dyoff = 0.5 * ydo + offconst * xdo
                if dxoff * dxoff + dyoff * dyoff < ddx_cc * ddx_cc + ddy_cc * ddy_cc:
                    tcx, tcy = a.x + dxoff, a.y + dyoff
            elif aodist < dadist:
                dxoff = 0.5 * xao + offconst * yao
                dyoff = 0.5 * yao - offconst * xao
                if dxoff * dxoff + dyoff * dyoff < ddx_cc * ddx_cc + ddy_cc * ddy_cc:
                    tcx, tcy = a.x + dxoff, a.y + dyoff
            else:
                dxoff = 0.5 * (c.x - b.x) - offconst * (c.y - b.y)
                dyoff = 0.5 * (c.y - b.y) + offconst * (c.x - b.x)
                if (dxoff * dxoff + dyoff * dyoff <
                    (ddx_cc - xdo) ** 2 + (ddy_cc - ydo) ** 2):
                    tcx, tcy = b.x + dxoff, b.y + dyoff
        else:
            tcx = (a.x + b.x + c.x) / 3.0
            tcy = (a.y + b.y + c.y) / 3.0

        # Pre-insertion encroachment check
        enc_seg = find_encroached(tcx, tcy)
        if enc_seg >= 0 and opts.no_steiner and mesh.segments[enc_seg].marker != 0:
            continue
        if enc_seg >= 0:
            seg_aff = []
            if split_segment(enc_seg, seg_aff, bad_idx):
                dbg_encroach += 1
                mr = tri_metric(mesh, bad_idx, min_ang_val, global_max_a, use_region_area, cos_min_ang)
                if mr[0] > 0:
                    heappush(pq, (-mr[0], bad_idx, mesh.triangles[bad_idx].generation))
                for ti2 in seg_aff:
                    if 0 <= ti2 < len(mesh.triangles) and mesh.triangles[ti2].v[0] >= 0:
                        for k in range(3):
                            nb = mesh.triangles[ti2].neighbors[k]
                            if nb >= 0:
                                mr2 = tri_metric(mesh, nb, min_ang_val, global_max_a,
                                                 use_region_area, cos_min_ang)
                                if mr2[0] > 0:
                                    heappush(pq, (-mr2[0], nb, mesh.triangles[nb].generation))
            continue

        # Insert
        np_pt = Point(tcx, tcy, len(mesh.vertices), 0)
        affected = []
        new_pt = insert_point_lawson(mesh, np_pt, bad_idx, constrained_edges, affected)
        if new_pt < 0:
            dbg_fail += 1
            continue

        steiner_count[0] += 1
        dbg_insert += 1

        # Use a flat sorted list instead of set to avoid hash overhead
        to_check = list(affected)
        for ti in affected:
            if 0 <= ti < len(mesh.triangles) and mesh.triangles[ti].v[0] >= 0:
                for j in range(3):
                    nb = mesh.triangles[ti].neighbors[j]
                    if nb >= 0:
                        to_check.append(nb)
        to_check.sort()
        # Deduplicate
        to_check = list(dict.fromkeys(to_check))
        for ti in to_check:
            if 0 <= ti < len(mesh.triangles):
                mr = tri_metric(mesh, ti, min_ang_val, global_max_a, use_region_area, cos_min_ang)
                if mr[0] > 0:
                    heappush(pq, (-mr[0], ti, mesh.triangles[ti].generation))

    if not opts.quiet:
        print(f"Quality refinement: {dbg_insert} Steiner points, "
              f"{dbg_encroach} segment splits, {dbg_skip} skipped, "
              f"{dbg_fail} failed insertions", file=sys.stderr)


# ============================================================
# Edge extraction
# ============================================================

def extract_edges(mesh):
    mesh.edges = []
    seen = set()
    for t in mesh.triangles:
        if t.v[0] < 0:
            continue
        for j in range(3):
            key = edge_key(t.v[j], t.v[(j + 1) % 3])
            if key not in seen:
                seen.add(key)
                mesh.edges.append((t.v[j], t.v[(j + 1) % 3]))


# ============================================================
# File I/O
# ============================================================

def skip_comments(lines, idx):
    while idx < len(lines) and (not lines[idx].strip() or lines[idx].strip().startswith('#')):
        idx += 1
    return idx


def read_tokens(lines, idx):
    idx = skip_comments(lines, idx)
    if idx >= len(lines):
        return [], idx
    return lines[idx].split(), idx + 1


def read_node_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 0
    tokens, idx = read_tokens(lines, idx)
    n, dim, n_attribs, n_markers = int(tokens[0]), int(tokens[1]), int(tokens[2]), int(tokens[3])
    pts = []
    for _ in range(n):
        tokens, idx = read_tokens(lines, idx)
        pt_idx = int(tokens[0])
        x, y = float(tokens[1]), float(tokens[2])
        attribs = [float(tokens[3 + a]) for a in range(n_attribs)]
        marker = int(tokens[3 + n_attribs]) if n_markers > 0 else 0
        lfs_col = 3 + n_attribs + (1 if n_markers > 0 else 0)
        lfs = float(tokens[lfs_col]) if len(tokens) > lfs_col else -1.0
        pts.append(Point(x, y, pt_idx, marker, attribs, lfs))
    return pts, n_attribs, n_markers


def read_poly_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 0
    tokens, idx = read_tokens(lines, idx)
    nv, dim, n_attribs, nmark = int(tokens[0]), int(tokens[1]), int(tokens[2]), int(tokens[3])
    pts = []
    first_vert_idx = float('inf')
    for _ in range(nv):
        tokens, idx = read_tokens(lines, idx)
        pt_idx = int(tokens[0])
        first_vert_idx = min(first_vert_idx, pt_idx)
        x, y = float(tokens[1]), float(tokens[2])
        attribs = [float(tokens[3 + a]) for a in range(n_attribs)]
        marker = int(tokens[3 + n_attribs]) if nmark > 0 else 0
        lfs_col = 3 + n_attribs + (1 if nmark > 0 else 0)
        lfs = float(tokens[lfs_col]) if len(tokens) > lfs_col else -1.0
        pts.append(Point(x, y, pt_idx, marker, attribs, lfs))
    base = 0 if first_vert_idx == 0 else 1

    tokens, idx = read_tokens(lines, idx)
    ns, smark = int(tokens[0]), int(tokens[1])
    segs = []
    for _ in range(ns):
        tokens, idx = read_tokens(lines, idx)
        v0, v1 = int(tokens[1]) - base, int(tokens[2]) - base
        marker = int(tokens[3]) if smark > 0 else 0
        lfs_col = 3 + (1 if smark > 0 else 0)
        lfs = float(tokens[lfs_col]) if len(tokens) > lfs_col else -1.0
        segs.append(Segment(v0, v1, marker, lfs))

    tokens, idx = read_tokens(lines, idx)
    nh = int(tokens[0])
    holes = []
    for _ in range(nh):
        tokens, idx = read_tokens(lines, idx)
        holes.append(Hole(float(tokens[1]), float(tokens[2])))

    regions = []
    if idx < len(lines):
        tokens, idx = read_tokens(lines, idx)
        if tokens:
            nr = int(tokens[0])
            for _ in range(nr):
                tokens, idx = read_tokens(lines, idx)
                rx, ry = float(tokens[1]), float(tokens[2])
                attrib = float(tokens[3])
                max_area = float(tokens[4]) if len(tokens) > 4 else -1.0
                regions.append(Region(rx, ry, attrib, max_area))

    # ---- tangle extensions: arcs and PBCs (optional, after regions) ----

    # Arc segments
    tokens, idx = read_tokens(lines, idx)
    if tokens:
        n_arcs = int(tokens[0])
        arc_has_markers = int(tokens[1]) if len(tokens) > 1 else 0
        for _ in range(n_arcs):
            tokens, idx = read_tokens(lines, idx)
            v0_raw, v1_raw = int(tokens[1]), int(tokens[2])
            v0a, v1a = v0_raw - base, v1_raw - base
            angle = float(tokens[3])
            max_seg_angle = float(tokens[4])
            col = 5
            arc_marker = 0
            if arc_has_markers and len(tokens) > col:
                arc_marker = int(tokens[col]); col += 1
            arc_lfs = float(tokens[col]) if len(tokens) > col else -1.0

            # Sanity checks
            angle = abs(angle)
            if angle < 1.0 or angle > 180.0:
                print(f"Warning: arc angle {angle}° clamped to [1,180]", file=sys.stderr)
                angle = max(1.0, min(180.0, angle))
            if max_seg_angle > 10.0 or max_seg_angle <= 0.01:
                max_seg_angle = 10.0

            # Discretize arc into chord segments
            ax, ay = pts[v0a].x, pts[v0a].y
            bx, by = pts[v1a].x, pts[v1a].y
            d = math.sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay))
            if d < 1e-30:
                segs.append(Segment(v0a, v1a, arc_marker, arc_lfs))
                continue

            tta = angle * math.pi / 180.0
            R = d / (2.0 * math.sin(tta / 2.0))
            tx, ty = (bx - ax) / d, (by - ay) / d
            h = math.sqrt(max(0.0, R * R - d * d / 4.0))
            cx = ax + (d / 2.0) * tx + h * (-ty)
            cy = ay + (d / 2.0) * ty + h * tx

            k = max(1, math.ceil(angle / max_seg_angle))

            if k == 1:
                segs.append(Segment(v0a, v1a, arc_marker, arc_lfs))
            else:
                dtheta = angle * math.pi / (180.0 * k)
                cos_d, sin_d = math.cos(dtheta), math.sin(dtheta)
                px, py = ax, ay
                prev_idx = v0a
                for j in range(k):
                    dx, dy = px - cx, py - cy
                    px = cx + dx * cos_d - dy * sin_d
                    py = cy + dx * sin_d + dy * cos_d
                    if j < k - 1:
                        new_idx = len(pts)
                        pts.append(Point(px, py, new_idx, arc_marker, [], arc_lfs))
                        segs.append(Segment(prev_idx, new_idx, arc_marker, arc_lfs))
                        prev_idx = new_idx
                    else:
                        segs.append(Segment(prev_idx, v1a, arc_marker, arc_lfs))

    # PBC definitions: pair two boundary markers
    tokens, idx = read_tokens(lines, idx)
    if tokens:
        n_pbcs = int(tokens[0])
        for _ in range(n_pbcs):
            tokens, idx = read_tokens(lines, idx)
            marker_a, marker_b = int(tokens[1]), int(tokens[2])
            pbc_type = int(tokens[3])
            for s in segs:
                if s.marker == marker_a or s.marker == marker_b:
                    s.pbc_type = pbc_type

    return pts, segs, holes, regions, n_attribs


def write_node_file(fn, pts, n_attribs, has_markers, opts):
    has_lfs = any(p.lfs > 0 for p in pts)
    with open(fn, 'w') as f:
        f.write(f"{len(pts)} 2 {n_attribs} {1 if has_markers else 0}\n")
        for i, pt in enumerate(pts):
            out_idx = i if opts.zero_indexed else i + 1
            line = f"{out_idx} {pt.x:.17g} {pt.y:.17g}"
            for a in pt.attribs:
                line += f" {a}"
            if has_markers:
                line += f" {pt.marker}"
            if has_lfs:
                line += f" {pt.lfs:.17g}"
            f.write(line + "\n")


def write_ele_file(fn, mesh, has_attrib, opts):
    with open(fn, 'w') as f:
        f.write(f"{len(mesh.triangles)} 3 {1 if has_attrib else 0}\n")
        for i, t in enumerate(mesh.triangles):
            out_idx = i if opts.zero_indexed else i + 1
            off = 0 if opts.zero_indexed else 1
            v0, v1, v2 = t.v[0] + off, t.v[1] + off, t.v[2] + off
            line = f"{out_idx} {v0} {v1} {v2}"
            if has_attrib:
                line += f" {t.region_attrib}"
            f.write(line + "\n")


def write_edge_file(fn, mesh, opts):
    with open(fn, 'w') as f:
        edge_count = {}
        for t in mesh.triangles:
            if t.v[0] < 0:
                continue
            for j in range(3):
                edge_count[edge_key(t.v[j], t.v[(j + 1) % 3])] = \
                    edge_count.get(edge_key(t.v[j], t.v[(j + 1) % 3]), 0) + 1
        seg_marker = {}
        for s in mesh.segments:
            seg_marker[edge_key(s.v0, s.v1)] = s.marker
        f.write(f"{len(mesh.edges)} 1\n")
        for i, (a, b) in enumerate(mesh.edges):
            out_idx = i if opts.zero_indexed else i + 1
            off = 0 if opts.zero_indexed else 1
            oa, ob = a + off, b + off
            key = edge_key(a, b)
            marker = 0
            is_boundary = edge_count.get(key, 0) == 1
            if key in seg_marker:
                marker = seg_marker[key]
                if marker == 0 and is_boundary:
                    marker = 1
            elif is_boundary:
                marker = 1
            f.write(f"{out_idx} {oa} {ob} {marker}\n")


def write_neigh_file(fn, mesh, opts):
    with open(fn, 'w') as f:
        f.write(f"{len(mesh.triangles)} 3\n")
        off = 0 if opts.zero_indexed else 1
        for i, t in enumerate(mesh.triangles):
            out_idx = i if opts.zero_indexed else i + 1

            def nb(n):
                return -1 if n < 0 else n + off

            f.write(f"{out_idx} {nb(t.neighbors[0])} {nb(t.neighbors[1])} {nb(t.neighbors[2])}\n")


def write_poly_file(fn, mesh, opts):
    has_lfs = any(p.lfs > 0 for p in mesh.vertices)
    has_seg_lfs = any(s.lfs > 0 for s in mesh.segments)
    with open(fn, 'w') as f:
        nv = len(mesh.vertices)
        f.write(f"{nv} 2 0 1\n")
        off = 0 if opts.zero_indexed else 1
        for i, v in enumerate(mesh.vertices):
            out_idx = i + off
            line = f"{out_idx} {v.x:.17g} {v.y:.17g} {v.marker}"
            if has_lfs:
                line += f" {v.lfs:.17g}"
            f.write(line + "\n")
        f.write(f"{len(mesh.segments)} 1\n")
        for i, s in enumerate(mesh.segments):
            out_idx = i + off
            line = f"{out_idx} {s.v0 + off} {s.v1 + off} {s.marker}"
            if has_seg_lfs:
                line += f" {s.lfs:.17g}"
            f.write(line + "\n")
        f.write("0\n")


def write_pbc_file(fn, mesh, opts):
    with open(fn, 'w') as f:
        f.write(f"{len(mesh.pbc_pairs)}\n")
        for i, p in enumerate(mesh.pbc_pairs):
            f.write(f"{i}\t{p.node_a}\t{p.node_b}\t{p.type}\n")



# ============================================================
# Argument parsing
# ============================================================

def parse_options(switch_str):
    opts = Options()
    i = 0
    while i < len(switch_str):
        c = switch_str[i]
        i += 1
        if c == 'p':
            opts.pslg = True
        elif c == 'P':
            opts.no_poly_out = True
        elif c == 'j':
            opts.jettison = True
        elif c == 'q':
            opts.quality = True
            num = ''
            while i < len(switch_str) and (switch_str[i].isdigit() or switch_str[i] == '.'):
                num += switch_str[i]
                i += 1
            if num:
                opts.min_angle = float(num)
        elif c == 'e':
            opts.edges = True
        elif c == 'A':
            opts.regions = True
        elif c == 'a':
            opts.area_limit = True
            num = ''
            while i < len(switch_str) and (switch_str[i].isdigit() or switch_str[i] == '.'):
                num += switch_str[i]
                i += 1
            if num:
                opts.max_area = float(num)
        elif c == 'z':
            opts.zero_indexed = True
            opts.first_index = 0
        elif c == 'Q':
            opts.quiet = True
        elif c == 'I':
            opts.suppress_iter = True
        elif c == 'O':
            opts.no_holes = True
        elif c == 'Y':
            opts.no_steiner = True
        elif c == 'n':
            opts.neighbors = True
        elif c == 'C':
            opts.clean_pslg = True
            num = ''
            while i < len(switch_str) and (switch_str[i].isdigit() or switch_str[i] == '.'):
                num += switch_str[i]
                i += 1
            if num:
                opts.clean_tol = float(num)
    return opts


# ============================================================
# Main
# ============================================================

def print_usage():
    print("""Usage: tangle.py [switches] inputfile

Switches:
  -p    Read a Planar Straight-Line Graph (.poly file) and triangulate it.
  -P    Suppress output of .poly file.
  -j    Jettison vertices not present in any triangle from .node output.
  -q    Quality mesh generation (min angle, default 20 deg, e.g. -q28.5).
  -e    Output .edge file.
  -A    Assign regional attributes from .poly regions.
  -a    Area constraints (optional global max, e.g. -a0.5).
  -z    Number items starting from zero.
  -Q    Quiet mode.
  -I    Suppress iteration numbers on output file names.
  -Y    No Steiner points on boundary segments.
  -n    Output .neigh file.

Input: .node or .poly files. Output: .1.node .1.ele [.1.edge] [.1.neigh] [.1.poly]""",
          file=sys.stderr)


def main():
    if len(sys.argv) < 2:
        print_usage()
        return 1

    switch_str = ''
    input_file = ''
    for arg in sys.argv[1:]:
        if arg.startswith('-'):
            switch_str += arg[1:]
        else:
            input_file = arg
    if not input_file:
        print("No input file specified.", file=sys.stderr)
        return 1

    opts = parse_options(switch_str)

    # Determine base and extension
    dot = input_file.rfind('.')
    if dot >= 0:
        base = input_file[:dot]
        ext = input_file[dot:]
    else:
        base = input_file
        ext = ''
    if not ext:
        ext = '.poly' if opts.pslg else '.node'
        input_file += ext

    mesh = Mesh()
    n_attribs = 0
    n_markers = 0

    if opts.pslg or ext == '.poly':
        opts.pslg = True
        pts, segs, holes, regions, n_attribs = read_poly_file(input_file)
        mesh.vertices = pts
        mesh.segments = segs
        mesh.holes = holes
        mesh.regions = regions

        # Build PBC twin map for .poly PBC definitions.
        # PBC in .poly uses marker pairs: segments with marker_a pair with marker_b.
        pbc_markers = set()
        for s in mesh.segments:
            if s.pbc_type >= 0:
                pbc_markers.add(s.marker)

        if pbc_markers:
            def build_chain(seg_indices):
                if not seg_indices:
                    return []
                node_segs = {}
                for si in seg_indices:
                    s = mesh.segments[si]
                    node_segs.setdefault(s.v0, []).append(si)
                    node_segs.setdefault(s.v1, []).append(si)
                start_node = -1
                for node, slist in node_segs.items():
                    if len(slist) == 1:
                        start_node = node; break
                if start_node < 0:
                    return []
                chain = [start_node]
                visited = set()
                cur_node = start_node
                while True:
                    found = False
                    for si in node_segs[cur_node]:
                        if si in visited:
                            continue
                        visited.add(si)
                        s = mesh.segments[si]
                        nxt = s.v1 if s.v0 == cur_node else s.v0
                        chain.append(nxt); cur_node = nxt; found = True; break
                    if not found:
                        break
                return chain

            # Group PBC markers by pbc_type into pairs
            type_to_markers = {}
            for m in pbc_markers:
                for s in mesh.segments:
                    if s.marker == m and s.pbc_type >= 0:
                        type_to_markers.setdefault(s.pbc_type, []).append(m)
                        break

            for pbc_type, markers in type_to_markers.items():
                for mi in range(0, len(markers) - 1, 2):
                    m_a, m_b = markers[mi], markers[mi + 1]
                    segs_a = [si for si, s in enumerate(mesh.segments)
                              if s.marker == m_a and s.pbc_type >= 0]
                    segs_b = [si for si, s in enumerate(mesh.segments)
                              if s.marker == m_b and s.pbc_type >= 0]
                    chain_a = build_chain(segs_a)
                    chain_b = build_chain(segs_b)
                    if not chain_a or not chain_b:
                        continue
                    if len(chain_a) != len(chain_b):
                        print(f"Warning: PBC markers {m_a} and {m_b} have different "
                              f"chain lengths ({len(chain_a)} vs {len(chain_b)})",
                              file=sys.stderr)
                        continue
                    # Orient chains consistently
                    def dist2(a, b):
                        dx = mesh.vertices[a].x - mesh.vertices[b].x
                        dy = mesh.vertices[a].y - mesh.vertices[b].y
                        return dx * dx + dy * dy
                    d_same = dist2(chain_a[0], chain_b[0]) + dist2(chain_a[-1], chain_b[-1])
                    d_rev = dist2(chain_a[0], chain_b[-1]) + dist2(chain_a[-1], chain_b[0])
                    if d_rev < d_same:
                        chain_b.reverse()
                    for i in range(len(chain_a)):
                        a, b = chain_a[i], chain_b[i]
                        if a == b:
                            continue
                        mesh.pbc_twin[a] = b; mesh.pbc_twin[b] = a
                        mesh.pbc_node_type[a] = pbc_type
                        mesh.pbc_node_type[b] = pbc_type
    else:
        pts, n_attribs, n_markers = read_node_file(input_file)
        mesh.vertices = pts

    # Shift coordinates to near the origin to maximize floating-point precision.
    shift_x, shift_y = 0.0, 0.0
    if mesh.vertices:
        min_x = min(p.x for p in mesh.vertices)
        max_x = max(p.x for p in mesh.vertices)
        min_y = min(p.y for p in mesh.vertices)
        max_y = max(p.y for p in mesh.vertices)
        shift_x = (min_x + max_x) * 0.5
        shift_y = (min_y + max_y) * 0.5
        for p in mesh.vertices:
            p.x -= shift_x; p.y -= shift_y
        for h in mesh.holes:
            h.x -= shift_x; h.y -= shift_y
        for r in mesh.regions:
            r.x -= shift_x; r.y -= shift_y

    if not opts.quiet:
        msg = f"Input: {len(mesh.vertices)} vertices"
        if mesh.segments:
            msg += f", {len(mesh.segments)} segments"
        if mesh.holes:
            msg += f", {len(mesh.holes)} holes"
        if mesh.regions:
            msg += f", {len(mesh.regions)} regions"
        print(msg, file=sys.stderr)

    t_start = time.perf_counter()
    t_prev = [t_start]

    def elapsed():
        now = time.perf_counter()
        ms = (now - t_prev[0]) * 1000
        t_prev[0] = now
        return ms

    # 0. Optional PSLG cleanup
    if opts.clean_pslg and opts.pslg:
        clean_pslg(mesh, opts.clean_tol, opts.quiet)
        if not opts.quiet:
            print(f"  PSLG cleanup: {elapsed():.1f} ms", file=sys.stderr)

    # 1. Delaunay triangulation
    build_delaunay(mesh)
    if not opts.quiet:
        print(f"Initial Delaunay: {len(mesh.triangles)} triangles ({elapsed():.1f} ms)",
              file=sys.stderr)

    # 2. Enforce PSLG segments (CDT)
    mesh.rebuild_adjacency()
    if opts.pslg:
        miss = enforce_constraints(mesh, opts.quiet)

        # Split long unenforced segments at their midpoint and rebuild.
        for split_round in range(10):
            if miss <= 0:
                break
            sum_seg_len = sum(
                math.sqrt((mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)*(mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)+
                           (mesh.vertices[s.v1].y-mesh.vertices[s.v0].y)*(mesh.vertices[s.v1].y-mesh.vertices[s.v0].y))
                for s in mesh.segments
            )
            avg_seg_len = sum_seg_len / max(1, len(mesh.segments))
            min_split_len = avg_seg_len * 2.0

            v2t_chk = [[] for _ in range(len(mesh.vertices))]
            for i, t in enumerate(mesh.triangles):
                if t.v[0] < 0:
                    continue
                for j in range(3):
                    v2t_chk[t.v[j]].append(i)

            to_split = []
            for i, s in enumerate(mesh.segments):
                if mesh.has_edge(s.v0, s.v1, v2t_chk):
                    continue
                seg_len = math.sqrt((mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)*(mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)+
                                     (mesh.vertices[s.v1].y-mesh.vertices[s.v0].y)*(mesh.vertices[s.v1].y-mesh.vertices[s.v0].y))
                if seg_len >= min_split_len:
                    to_split.append(i)
            if not to_split:
                break

            for si in reversed(to_split):
                seg = mesh.segments[si]
                mx = (mesh.vertices[seg.v0].x + mesh.vertices[seg.v1].x) / 2
                my = (mesh.vertices[seg.v0].y + mesh.vertices[seg.v1].y) / 2
                bm = max(mesh.vertices[seg.v0].marker, mesh.vertices[seg.v1].marker)
                mid_idx = len(mesh.vertices)
                mesh.vertices.append(Point(mx, my, mid_idx, bm, []))
                v0, v1, mk = seg.v0, seg.v1, seg.marker
                mesh.segments[si] = Segment(v0, mid_idx, mk)
                mesh.segments.insert(si + 1, Segment(mid_idx, v1, mk))

            if not opts.quiet:
                print(f"  Split {len(to_split)} long unenforced segments, rebuilding...",
                      file=sys.stderr)
            mesh.triangles.clear()
            build_delaunay(mesh)
            mesh.rebuild_adjacency()
            miss = enforce_constraints(mesh, opts.quiet)

    if not opts.quiet:
        print(f"  CDT: {elapsed():.1f} ms", file=sys.stderr)

    # 3. Remove holes
    if opts.pslg:
        remove_holes(mesh, opts)
    if not opts.quiet:
        print(f"  Holes: {elapsed():.1f} ms", file=sys.stderr)

    # 4. Assign regions
    if opts.regions:
        assign_regions(mesh, opts)
    if not opts.quiet:
        print(f"  Regions: {elapsed():.1f} ms", file=sys.stderr)

    n_input_verts = len(mesh.vertices)

    # 5. Quality refinement (new triangles inherit regions from containing triangle)
    if opts.quality or opts.area_limit:
        refine_quality(mesh, opts, n_input_verts)
    if not opts.quiet:
        print(f"  Refinement: {elapsed():.1f} ms", file=sys.stderr)

    # 7. Final cleanup: compact dead triangles, rebuild adjacency, extract edges
    mesh.triangles = [t for t in mesh.triangles if t.v[0] >= 0]
    mesh.rebuild_adjacency()
    extract_edges(mesh)
    if not opts.quiet:
        print(f"  Cleanup+edges: {elapsed():.1f} ms", file=sys.stderr)

    # Convert PBC twin map to output pairs
    seen_pbc = set()
    for a, b in mesh.pbc_twin.items():
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key not in seen_pbc:
            seen_pbc.add(key)
            pbc_type = mesh.pbc_node_type.get(a, 0)
            mesh.pbc_pairs.append(PBCNodePair(a, b, pbc_type))

    # Jettison
    if opts.jettison:
        used = [False] * len(mesh.vertices)
        for t in mesh.triangles:
            used[t.v[0]] = used[t.v[1]] = used[t.v[2]] = True
        vremap = [-1] * len(mesh.vertices)
        new_pts = []
        for i in range(len(mesh.vertices)):
            if used[i]:
                vremap[i] = len(new_pts)
                new_pts.append(mesh.vertices[i])
        mesh.vertices = new_pts
        for t in mesh.triangles:
            t.v = [vremap[t.v[0]], vremap[t.v[1]], vremap[t.v[2]]]
        mesh.edges = [(vremap[a], vremap[b]) for a, b in mesh.edges]
        for s in mesh.segments:
            s.v0 = vremap[s.v0]
            s.v1 = vremap[s.v1]
        for p in mesh.pbc_pairs:
            p.node_a = vremap[p.node_a]
            p.node_b = vremap[p.node_b]

    if not opts.quiet:
        print(f"Output: {len(mesh.vertices)} vertices, {len(mesh.triangles)} triangles",
              file=sys.stderr)

    # Shift coordinates back to original position for output
    for p in mesh.vertices:
        p.x += shift_x; p.y += shift_y

    out_base = base + ("" if opts.suppress_iter else ".1")
    write_node_file(out_base + ".node", mesh.vertices, n_attribs,
                    n_markers > 0 or opts.pslg, opts)
    write_ele_file(out_base + ".ele", mesh, opts.regions, opts)
    if opts.edges:
        write_edge_file(out_base + ".edge", mesh, opts)
    if opts.neighbors:
        write_neigh_file(out_base + ".neigh", mesh, opts)
    if opts.pslg and not opts.no_poly_out:
        write_poly_file(out_base + ".poly", mesh, opts)
    if mesh.pbc_pairs:
        write_pbc_file(out_base + ".pbc", mesh, opts)

    if not opts.quiet:
        msg = f"Wrote {out_base}.node, {out_base}.ele"
        if opts.edges:
            msg += f", {out_base}.edge"
        if opts.neighbors:
            msg += f", {out_base}.neigh"
        if opts.pslg and not opts.no_poly_out:
            msg += f", {out_base}.poly"
        if mesh.pbc_pairs:
            msg += f", {out_base}.pbc"
        print(msg, file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main())
