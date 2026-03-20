# tangle

**Version 0.2** — 17 Mar 2026

A 2D constrained Delaunay triangulation and quality mesh generation tool,
implemented in C++17. File format compatible with Shewchuk's
[Triangle](https://www.cs.cmu.edu/~quake/triangle.html) with extensions
for arc segments, local feature size constraints, and periodic boundary
conditions.

A Python implementation (`tangle.py`) is also included with identical
functionality.

## Authors

- **David Meeker** — author
- **Claude Code (Opus 4.6)** — development assistance

## Why tangle?

Tangle is backward-compatible with Triangle's `.node`, `.poly`, `.ele`, `.edge`, and `.neigh` file formats for compatibility with existing workflows. Tangle also provides additional capabilities that are practically important to formulation of meshes for the solution of finite element problems (particularly in the area of electric machines):

- **Local feature size (LFS) constraints** — per-vertex and per-segment
  mesh density control, so that refinement can be concentrated where it
  matters without over-meshing the rest of the domain.
- **Periodic boundary conditions (PBCs)** — paired boundary chains with
  synchronized segment splitting during refinement, producing matching
  node-to-node correspondence for periodic or anti-periodic coupling.
- **Arc segments** — circular arcs defined by two endpoints and a
  subtended angle, automatically discretized into chord segments.

These extensions are described in [`polyformat.txt`](https://www.femm.info/wiki/PolyFormat),
which documents the full extended `.poly` format.

## Build

```bash
g++ -O3 -std=c++17 -static -o tangle tangle.cpp float256.cpp -lm
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
| `-C`   | PSLG cleanup: merge near-duplicate nodes, split intersecting segments, remove duplicates. Optional tolerance (e.g. `-C0.001`); default is bounding-box diagonal × 1e-6 |
| `-Y`   | No Steiner points on boundary segments |
| `-Q`   | Quiet mode |
| `-I`   | Suppress iteration numbers on output filenames |
| `-j`   | Jettison unused vertices from output |
| `-P`   | Suppress .poly file output |

**Input:** `.node` (point set) or `.poly` (planar straight-line graph) files,
using the same format as Triangle.

**Output:** `<base>.1.node`, `<base>.1.ele`, and optionally `.edge`, `.neigh`,
`.poly`, `.pbc`.

### Example

```bash
# Generate a quality mesh from a PSLG with 30-degree minimum angle
tangle -pq30 myExample

# Same, with area constraint and edge output
tangle -pq30 -a0.1 -e myExample
```

## Algorithms

Tangle's pipeline has four main stages, followed by extension processing
for LFS, arcs, and PBCs.

### 1. Delaunay Triangulation — Bowyer-Watson Algorithm

The initial triangulation uses the [Bowyer-Watson algorithm](https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm) ([Bowyer 1981](https://adrianbowyer.com/Publications/computing-dirichlet-tessellations.pdf), [Watson 1981](https://doi.org/10.1093/comjnl/24.2.167)), an incremental insertion method. Points are inserted one at a time into a super-triangle that encloses the entire point set. For each new point, the algorithm finds all existing triangles whose circumcircle contains the point (the "cavity"), removes them, and re-triangulates the resulting polygonal hole with triangles connecting the new point to each edge of the cavity boundary. Points are sorted by x-coordinate before insertion so that a walking point-location strategy (starting from the last inserted triangle) is efficient.

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

Holes are removed by Breadth-First Search (BFS) flood-fill from seed points
specified in the .poly file, stopping at constrained edges. Region attributes
and per-region area constraints are assigned by a similar flood-fill from
region seed points.

### 4. Quality Refinement — Ruppert's Algorithm with Üngör Off-Centers

Bad triangles (minimum angle below the threshold, or area above the constraint) are refined using [Ruppert's Delaunay refinement algorithm](https://en.wikipedia.org/wiki/Ruppert%27s_algorithm) ([Ruppert 1995](https://www.cis.upenn.edu/~cis6100/ruppert.pdf)) with Üngör off-center point placement ([Üngör 2004](https://www.cise.ufl.edu/~ungor/papers/offCenter.pdf)). Instead of inserting at the circumcenter, the off-center method places the new vertex closer to the shortest edge of the bad triangle, producing better-shaped elements and faster convergence.

Before inserting a new vertex, the algorithm checks whether it would
encroach on any constrained segment (i.e., fall inside the segment's
diametral circle). If so, the segment is split at its midpoint instead.
This guarantees termination: segment splitting resolves encroachment, and
the off-center placement ensures that new vertices do not create
progressively worse triangles.

New vertices are inserted into the existing triangulation using [Lawson's incremental insertion algorithm](https://en.wikipedia.org/wiki/Delaunay_triangulation#Flip_algorithms) ([Lawson 1977](https://ntrs.nasa.gov/api/citations/19770025881/downloads/19770025881.pdf)), which splits the containing triangle and restores the Delaunay property through local edge flips.

A priority queue processes the worst triangles first. An adaptive Steiner
point limit, estimated from the constrained Delaunay triangulation using a
local feature size integral (Ruppert 1995), prevents non-convergence on
pathological inputs.

### Local Feature Size (LFS) Constraints

LFS values can be specified per-vertex and per-segment in the `.poly` file
(see [`polyformat.txt`](https://www.femm.info/wiki/PolyFormat) for syntax).  During quality
refinement, a triangle is considered "bad" not only if its minimum angle or
area violates the global thresholds, but also if any of its edges exceeds the
LFS value inherited from its endpoints.  This provides spatially varying
mesh density control: regions near vertices or segments with small LFS values
get fine elements, while the rest of the domain stays coarse.

### Arc Segments

Arc segments define circular arcs between two existing vertices by specifying
the subtended angle and a maximum discretization angle.  During input
processing, each arc is replaced by a chain of straight constrained segments
(chord segments) connecting points computed along the arc.  Interior vertices
created by the discretization are added to the point set and participate
normally in the Delaunay triangulation and all subsequent stages.

The arc center is determined from the chord and subtended angle, with the
center placed to the left of the chord as traversed from endpoint 1 to
endpoint 2.  The number of chord segments is `ceil(angle / max_seg_angle)`.

### Periodic Boundary Conditions (PBCs)

PBC definitions pair two boundary chains identified by their segment
boundary markers.  During quality refinement, PBC-paired segments are split
synchronously: when the refinement algorithm splits a segment on one chain
(e.g. because a new Steiner point encroaches on it), it also splits the
partner segment on the paired chain at the corresponding parametric position.
This maintains a one-to-one node correspondence between the two boundaries
throughout refinement.

After meshing, the final node pairs are written to a `.pbc` output file.
This is useful for assembling finite element systems with periodic or
anti-periodic coupling between boundaries — for example, exploiting
rotational symmetry in electric machines, or using the Kelvin transformation
to map unbounded exterior domains onto a finite mesh.

### Geometric Robustness — Adaptive Exact Arithmetic

The geometric predicates `orient2d` (point-line orientation) and `inCircle` (point-in-circumcircle) use a floating-point filter with exact arithmetic fallback. Each predicate first evaluates a fast double-precision (~53-bit) computation, then checks the result against [a forward error bound derived from the operand magnitudes](https://www.femm.info/wiki/PredicateBounds) ([Higham 2002](https://www.google.com/books/edition/Accuracy_and_Stability_of_Numerical_Algo/epilvM5MMxwC)). When the result is too close to zero relative to the bound to determine the sign reliably, the computation falls back to extended-precision arithmetic for an exact answer.

Two custom floating-point types provide the extended precision, both implemented
as one IEEE double (carrying the sign, exponent, and top 52 mantissa bits)
plus additional unsigned 64-bit integers (`uint64_t`) that extend the mantissa:

- **`float128`** — double + 1×uint64_t, giving ~117 bits of mantissa
  precision. Used for `orient2d`, where the 2×2 determinant structure
  requires at most ~106 bits for an exact sign determination.
- **`float256`** — double + 3×uint64_t, giving ~245 bits of mantissa
  precision. Used for `inCircle`, where the 3×3 determinant with squared
  terms requires up to ~212 bits.

This two-tier approach (double → float128 or float256) avoids the robustness
failures that plague naive floating-point implementations while keeping the
common case fast. Over 99% of predicate calls are resolved at double
precision, and the few that fall through use the narrowest exact type
sufficient for the predicate.

## File Formats

Tangle uses the same file formats as Triangle for `.node`, `.poly`, `.ele`,
`.edge`, and `.neigh` files. See the
[Triangle documentation](https://www.cs.cmu.edu/~quake/triangle.html) for
the base format, and [`polyformat.txt`](https://www.femm.info/wiki/PolyFormat) for tangle's
extensions (LFS columns, arc segments, PBC definitions, and `.pbc` output).

## License

Tangle is an independent implementation distributed under the MIT License. See [LICENSE](LICENSE) for details.
