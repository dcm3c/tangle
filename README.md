# tangle

**Version 0.1** — 13 Mar 2026

A 2D constrained Delaunay triangulation and quality mesh generation tool,
implemented in readable C++17 (~1100 lines). File-format compatible with
Shewchuk's [Triangle](https://www.cs.cmu.edu/~quake/triangle.html), so it
works as a drop-in replacement in existing workflows, e.g.
[FEMM](https://www.femm.info) 

A Python implementation (`tangle.py`) is also included with identical
functionality.

## Authors

- **David Meeker** — author
- **Claude Code (Opus 4.6)** — development assistance

## Why tangle?

Tangle is an independent implementation that achieves file-format
interoperability with Triangle through original code — it contains no
Triangle source code. It is MIT-licensed and can be freely used, modified,
and embedded in any project. It is also deliberately kept small and readable
so that it can be understood, debugged, and extended without difficulty.

## Build

```bash
g++ -O2 -std=c++17 -static -o tangle tangle.cpp float256.cpp -lm
```

No external dependencies. Requires a C++17 compiler (GCC, Clang, or MSVC).

Pre-built Windows binaries are available at
[https://www.femm.info/wiki/Tangle](https://www.femm.info/wiki/Tangle).

## Usage

```bash
tangle [switches] inputfile
```

| Switch | Description |
|--------|-------------|
| `-p`   | Input is a Planar Straight-Line Graph, or PSLG (.poly file) |
| `-q`   | Quality mesh generation (default 20 deg minimum angle) |
| `-q30` | Quality mesh with 30 degree minimum angle |
| `-a0.5`| Maximum triangle area constraint |
| `-A`   | Assign region attributes from .poly region list |
| `-e`   | Output edges (.edge file) |
| `-n`   | Output triangle neighbors (.neigh file) |
| `-z`   | Zero-indexed output (default is 1-indexed) |
| `-Y`   | No Steiner points on boundary segments |
| `-Q`   | Quiet mode |
| `-I`   | Suppress iteration numbers on output filenames |
| `-j`   | Jettison unused vertices from output |
| `-P`   | Suppress .poly file output |

**Input:** `.node` (point set) or `.poly` (planar straight-line graph) files,
using the same format as Triangle.

**Output:** `<base>.1.node`, `<base>.1.ele`, and optionally `.edge`, `.neigh`,
`.poly`.

### Example

```bash
# Generate a quality mesh from a PSLG with 30-degree minimum angle
tangle -pq30 myExample

# Same, with area constraint and edge output
tangle -pq30 -a0.1 -e myExample
```

## Algorithms

Tangle's pipeline has four main stages:

### 1. Delaunay Triangulation — Bowyer-Watson Algorithm

The initial triangulation uses the
[Bowyer-Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm)
(Bowyer 1981, Watson 1981), an incremental insertion method. Points are
inserted one at a time into a super-triangle that encloses the entire point
set. For each new point, the algorithm finds all existing triangles whose
circumcircle contains the point (the "cavity"), removes them, and
re-triangulates the resulting polygonal hole with triangles connecting the new
point to each edge of the cavity boundary. Points are sorted by x-coordinate
before insertion so that a walking point-location strategy (starting from the
last inserted triangle) is efficient.

### 2. Constrained Delaunay Triangulation (CDT)

After the initial Delaunay triangulation, PSLG segments that are not already
present as edges are enforced to produce a
[constrained Delaunay triangulation](https://en.wikipedia.org/wiki/Constrained_Delaunay_triangulation).
The algorithm identifies triangles crossed by each missing segment, removes
them, and re-triangulates the two polygonal regions on either side using a
sweep approach. Segments that cannot be enforced in one pass are split at
their midpoint and the process repeats. Constrained edges are then respected
as barriers during all subsequent operations.

### 3. Hole and Region Processing

Holes are removed by Breadth-First Search (BFS) flood-fill from seed points specified in the .poly
file, stopping at constrained edges. Region attributes and per-region area
constraints are assigned by a similar flood-fill from region seed points.

### 4. Quality Refinement — Ruppert's Algorithm with Üngör Off-Centers

Bad triangles (minimum angle below the threshold, or area above the
constraint) are refined using
[Ruppert's Delaunay refinement algorithm](https://en.wikipedia.org/wiki/Ruppert%27s_algorithm)
(Ruppert 1995) with Üngör off-center point placement (Üngör 2004). Instead
of inserting at the circumcenter, the off-center method places the new vertex
closer to the shortest edge of the bad triangle, producing better-shaped
elements and faster convergence.

Before inserting a new vertex, the algorithm checks whether it would
encroach on any constrained segment (i.e., fall inside the segment's
diametral circle). If so, the segment is split at its midpoint instead.
This guarantees termination: segment splitting resolves encroachment, and
the off-center placement ensures that new vertices do not create
progressively worse triangles.

New vertices are inserted into the existing triangulation using
[Lawson's incremental insertion algorithm](https://en.wikipedia.org/wiki/Delaunay_triangulation#Flip_algorithms)
(Lawson 1977), which splits the containing triangle and restores the Delaunay
property through local edge flips.

A priority queue processes the worst triangles first. An adaptive Steiner
point limit, estimated from the constrained Delaunay triangulation using a
local feature size integral (Ruppert 1995), prevents non-convergence on
pathological inputs.

### Geometric Robustness — Adaptive Exact Arithmetic

The geometric predicates `orient2d` (point-line orientation) and `inCircle`
(point-in-circumcircle) use an adaptive precision approach inspired by
Shewchuk's robust predicates (Shewchuk 1997). A fast double-precision test
with an error bound handles the common case; when the result is ambiguous, the
computation falls back to exact arithmetic via a custom quad-double type
(`float256`). This avoids the robustness failures that plague naive
floating-point implementations while keeping the common case fast. Coordinates
are shifted to the bounding-box centroid before processing to maximize
floating-point precision.

## File Formats

Tangle uses the same file formats as Triangle. See the
[Triangle documentation](https://www.cs.cmu.edu/~quake/triangle.html) for
format details.

## License

MIT License. See [LICENSE](LICENSE) for details.
