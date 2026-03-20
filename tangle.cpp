////////////////////////////////////////////////////
// tangle
// A 2D Delaunay triangulation tool compatible with Shewchuk's Triangle format.
//
// Author: David Meeker
// Generated with the assistance of Claude Code Opus 4.6
//
// Version 0.2
// 17 Mar 2026
//
// Supports: -p -P -j -q -e -A -a -z -Q -I -Y options
//
// Build: g++ -O3 -std=c++17 -o tangle tangle.cpp float256.cpp -lm
//
// File format compatibility:
//   .node  - vertex files
//   .poly  - PSLG (Planar Straight-Line Graph) files
//   .ele   - triangle output files
//   .edge  - edge output files
//   .neigh - neighbor output files
//
// MIT License
//
// Copyright (c) 2026 David Meeker
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <limits>
#include <cstring>
#include "float256.h"
#include <numeric>
#include <queue>


struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<long long>()(((long long)p.first<<32)|((unsigned int)p.second));
    }
};
using EdgeSet = std::unordered_set<std::pair<int,int>, PairHash>;

// ============================================================
// Constants and Utilities
// ============================================================

// Base relative tolerance.  EPS is scaled to the bounding-box diagonal in
// buildDelaunay(); EPS_SCALE controls that ratio and is also used directly
// for dimensionless comparisons.
static const double EPS_SCALE = 1e-10;
static double EPS = EPS_SCALE;

struct Options {
    bool pslg          = false;
    bool no_poly_out   = false;
    bool jettison      = false;
    bool quality       = false;
    double min_angle   = 20.0;
    bool edges         = false;
    bool regions       = false;
    bool area_limit    = false;
    double max_area    = -1.0;
    bool zero_indexed  = false;
    bool quiet         = false;
    bool suppress_iter = false;
    bool no_holes      = false;
    bool no_steiner    = false;
    bool neighbors     = false;
    bool clean_pslg    = false;
    double clean_tol   = -1.0;  // -1 = auto (bboxDiag * 1e-6)
    int  first_index   = 1;
};

// ============================================================
// Geometry primitives
// ============================================================

struct Point {
    double x, y;
    int    id;
    int    marker;
    std::vector<double> attribs;
    double lfs = -1.0;  // local feature size constraint (-1 = none)
};

struct Segment {
    int v0, v1;
    int marker;
    double lfs = -1.0;   // per-segment LFS constraint (-1 = none)
    int pbc_partner = -1; // index of paired PBC segment (-1 = none)
    int pbc_type = -1;    // 0=periodic, 1=anti-periodic (-1 = none)
};

struct Hole { double x, y; };

struct Region {
    double x, y;
    double attrib;
    double max_area; // -1 means no constraint
};

struct Triangle {
    std::array<int,3> v;
    std::array<int,3> neighbors; // -1 = boundary
    double region_attrib;
    double region_max_area;
    int    marker;
    int    generation = 0;
};

// Lightweight orient2d for triangle walks — no exact fallback needed because
// the walk is intrinsically robust to occasional sign errors (it self-corrects
// or falls through to the linear scan).  Also avoids constructing a Point with
// its heap-allocated attribs vector.
inline double orient2d_xy(double ax, double ay, double bx, double by, double cx, double cy){
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

inline double orient2d(const Point& a, const Point& b, const Point& c){
    double Ax = b.x - a.x, Ay = b.y - a.y;
    double Bx = c.x - a.x, By = c.y - a.y;
    double det = Ax * By - Bx * Ay;

    // Error bound: det terms are O(M²), check if result is well above fp noise
    double M = std::max({std::abs(Ax), std::abs(Ay), std::abs(Bx), std::abs(By)});
    if(std::abs(det) > 6.662e-16 * M * M) return det;

    // Fall back to exact float128 arithmetic (~117 bits, sufficient for orient2d)
    float128 fAx = float128(b.x) - float128(a.x), fAy = float128(b.y) - float128(a.y);
    float128 fBx = float128(c.x) - float128(a.x), fBy = float128(c.y) - float128(a.y);
    float128 fdet = fAx * fBy - fBx * fAy;
    return fdet.hi;
}

inline double inCircle(const Point& a, const Point& b, const Point& c, const Point& d){
    // 3x3 formulation: subtract d to eliminate one row
    double Ax = a.x - d.x, Ay = a.y - d.y;
    double Bx = b.x - d.x, By = b.y - d.y;
    double Cx = c.x - d.x, Cy = c.y - d.y;
    double Bs = Bx*Bx + By*By;
    double Cs = Cx*Cx + Cy*Cy;

    double det = Ax*(By*Cs - Cy*Bs)
               - Ay*(Bx*Cs - Cx*Bs)
               + (Ax*Ax + Ay*Ay)*(Bx*Cy - Cx*By);

    // Error bound: det terms are O(M^4), check if result is well above fp noise
    double M = std::max({std::abs(Ax), std::abs(Ay), std::abs(Bx), std::abs(By),
                         std::abs(Cx), std::abs(Cy)});
    if(std::abs(det) > 1.333e-14 * M * M * M * M) return det;

    // Fall back to exact float256 arithmetic (standard 3x3 formulation)
    float256 fax = float256(a.x) - float256(d.x), fay = float256(a.y) - float256(d.y);
    float256 fbx = float256(b.x) - float256(d.x), fby = float256(b.y) - float256(d.y);
    float256 fcx = float256(c.x) - float256(d.x), fcy = float256(c.y) - float256(d.y);
    float256 fbs = fbx*fbx + fby*fby;
    float256 fcs = fcx*fcx + fcy*fcy;
    float256 fdet = fax*(fby*fcs - fcy*fbs)
                  - fay*(fbx*fcs - fcx*fbs)
                  + (fax*fax + fay*fay)*(fbx*fcy - fcx*fby);
    return fdet.hi;
}

double triArea(const Point& a, const Point& b, const Point& c){
    return 0.5*std::abs(orient2d(a,b,c));
}

// Angle at vertex a in triangle (a, b, c), in degrees
double vertexAngle(const Point& a, const Point& b, const Point& c){
    double abx=b.x-a.x, aby=b.y-a.y;
    double acx=c.x-a.x, acy=c.y-a.y;
    double dot=abx*acx+aby*acy;
    double cross=abx*acy-aby*acx;
    return std::atan2(std::abs(cross),dot)*180.0/M_PI;
}

double minAngle(const Point& a, const Point& b, const Point& c){
    double la = std::sqrt((b.x-c.x)*(b.x-c.x)+(b.y-c.y)*(b.y-c.y));
    double lb = std::sqrt((a.x-c.x)*(a.x-c.x)+(a.y-c.y)*(a.y-c.y));
    double lc = std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
    if(la<EPS||lb<EPS||lc<EPS) return 0.0;
    double cosA = std::clamp((lb*lb+lc*lc-la*la)/(2*lb*lc),-1.0,1.0);
    double cosB = std::clamp((la*la+lc*lc-lb*lb)/(2*la*lc),-1.0,1.0);
    double cosC = std::clamp((la*la+lb*lb-lc*lc)/(2*la*lb),-1.0,1.0);
    return std::min({std::acos(cosA),std::acos(cosB),std::acos(cosC)})*180.0/M_PI;
}

void circumcenter(const Point& a, const Point& b, const Point& c,
                  double& ccx, double& ccy){
    double ax=b.x-a.x, ay=b.y-a.y;
    double bx=c.x-a.x, by=c.y-a.y;
    double D = 2.0*(ax*by - ay*bx);
    if(std::abs(D)<EPS){ ccx=(a.x+b.x+c.x)/3.0; ccy=(a.y+b.y+c.y)/3.0; return; }
    double ux = (by*(ax*ax+ay*ay) - ay*(bx*bx+by*by))/D;
    double uy = (ax*(bx*bx+by*by) - bx*(ax*ax+ay*ay))/D;
    ccx = a.x+ux; ccy = a.y+uy;
}

bool pointOnSegment(const Point& p, const Point& a, const Point& b){
    if(orient2d(a,b,p)!=0) return false;
    double dx=b.x-a.x, dy=b.y-a.y;
    double t;
    if(std::abs(dx)>EPS) t=(p.x-a.x)/dx;
    else if(std::abs(dy)>EPS) t=(p.y-a.y)/dy;
    else return false;
    return t>=-EPS_SCALE && t<=1+EPS_SCALE;  // dimensionless: fixed tolerance
}

static inline std::pair<int,int> edgeKey(int a, int b){
    return {std::min(a,b), std::max(a,b)};
}

// ============================================================
// Mesh structure with fast adjacency
// ============================================================

struct PBCNodePair {
    int node_a, node_b;
    int type; // 0=periodic, 1=anti-periodic
};

struct Mesh {
    std::vector<Point>    vertices;
    std::vector<Triangle> triangles;
    std::vector<Segment>  segments;
    std::vector<Hole>     holes;
    std::vector<Region>   regions;
    std::vector<std::pair<int,int>> edges;
    std::vector<PBCNodePair> pbc_pairs; // filled after refinement
    std::map<int,int> pbc_twin;          // node↔twin mapping for PBC boundaries
    std::map<int,int> pbc_node_type;     // node→PBC type (0=periodic, 1=anti-periodic)

    void rebuildAdjacency(){
        int nt=(int)triangles.size();
        std::unordered_map<long long, int> edgeToTriEdge;
        for(auto& t : triangles) t.neighbors = {-1,-1,-1};
        for(int i=0; i<nt; i++){
            if(triangles[i].v[0]<0) continue;
            for(int j=0;j<3;j++){
                int a=triangles[i].v[j], b=triangles[i].v[(j+1)%3];
                int mn=std::min(a,b), mx=std::max(a,b);
                long long key = ((long long)mn << 32) | (unsigned int)mx;
                auto it = edgeToTriEdge.find(key);
                if(it != edgeToTriEdge.end()){
                    int packed = it->second;
                    int ot = packed >> 2, oe = packed & 3;
                    triangles[i].neighbors[j] = ot;
                    triangles[ot].neighbors[oe] = i;
                    edgeToTriEdge.erase(it);
                } else {
                    edgeToTriEdge[key] = (i << 2) | j;
                }
            }
        }
    }

    int locateTriangle(double px, double py, int hint = 0) const {
        if(triangles.empty()) return -1;
        if(hint < 0 || hint >= (int)triangles.size()) hint = 0;
        int cur = hint;
        for(int step=0; step<(int)triangles.size()*2; step++){
            const auto& t = triangles[cur];
            if(t.v[0]<0){cur=(cur+1)%(int)triangles.size();continue;}
            const auto &v0=vertices[t.v[0]], &v1=vertices[t.v[1]], &v2=vertices[t.v[2]];
            double o0 = orient2d_xy(v0.x,v0.y, v1.x,v1.y, px,py);
            double o1 = orient2d_xy(v1.x,v1.y, v2.x,v2.y, px,py);
            double o2 = orient2d_xy(v2.x,v2.y, v0.x,v0.y, px,py);
            if(o0>=-EPS && o1>=-EPS && o2>=-EPS) return cur;
            int worst=-1; double wv=0;
            if(o0<wv){wv=o0;worst=0;} if(o1<wv){wv=o1;worst=1;} if(o2<wv){wv=o2;worst=2;}
            if(worst>=0 && t.neighbors[worst]>=0) cur=t.neighbors[worst];
            else break;
        }
        // Fallback: linear scan
        for(int i=0;i<(int)triangles.size();i++){
            const auto& t=triangles[i];
            if(t.v[0]<0) continue;
            const auto &v0=vertices[t.v[0]], &v1=vertices[t.v[1]], &v2=vertices[t.v[2]];
            if(orient2d_xy(v0.x,v0.y, v1.x,v1.y, px,py)>=-EPS &&
               orient2d_xy(v1.x,v1.y, v2.x,v2.y, px,py)>=-EPS &&
               orient2d_xy(v2.x,v2.y, v0.x,v0.y, px,py)>=-EPS)
                return i;
        }
        return -1;
    }

    // Check if edge (u,v) exists using v2t map
    bool hasEdge(int u, int v, const std::vector<std::vector<int>>& v2t) const {
        for(int ti : v2t[u]){
            const auto& t = triangles[ti];
            if(t.v[0]<0) continue;
            for(int j=0;j<3;j++){
                if((t.v[j]==u && t.v[(j+1)%3]==v)||(t.v[j]==v && t.v[(j+1)%3]==u))
                    return true;
            }
        }
        return false;
    }
};

// ============================================================
// Bowyer-Watson Delaunay triangulation
// ============================================================

void buildDelaunay(Mesh& mesh){
    auto& pts = mesh.vertices;
    int n = (int)pts.size();

    double minX=pts[0].x, maxX=pts[0].x, minY=pts[0].y, maxY=pts[0].y;
    for(auto& p : pts){
        minX=std::min(minX,p.x); maxX=std::max(maxX,p.x);
        minY=std::min(minY,p.y); maxY=std::max(maxY,p.y);
    }
    double dx=maxX-minX, dy=maxY-minY, d=std::max(dx,dy)*10.0;

    // Scale EPS to the bounding-box diagonal so tolerances stay meaningful
    // regardless of coordinate magnitude.
    double bboxDiag = std::sqrt(dx*dx + dy*dy);
    EPS = bboxDiag * EPS_SCALE;
    pts.push_back({minX-d, minY-3*d, n,   0, {}});
    pts.push_back({minX+3*d, minY-d, n+1, 0, {}});
    pts.push_back({minX-d, minY+3*d, n+2, 0, {}});

    auto& tris = mesh.triangles;
    tris.clear();
    tris.push_back(Triangle{{n,n+1,n+2},{-1,-1,-1},0,-1,0});

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b){
        return pts[a].x < pts[b].x || (pts[a].x==pts[b].x && pts[a].y<pts[b].y);
    });

    int lastTri = 0;
    std::vector<int> visitEpoch;
    int epoch = 0;

    for(int ii=0; ii<n; ii++){
        int pidx = order[ii];
        const Point& p = pts[pidx];

        // Locate containing triangle by walking
        int startTri = -1;
        {
            int cur = std::min(lastTri, (int)tris.size()-1);
            for(int step=0; step<(int)tris.size()*2; step++){
                const auto& t = tris[cur];
                if(t.v[0]<0){cur=0;continue;}
                double o0=orient2d(pts[t.v[0]],pts[t.v[1]],p);
                double o1=orient2d(pts[t.v[1]],pts[t.v[2]],p);
                double o2=orient2d(pts[t.v[2]],pts[t.v[0]],p);
                if(o0>=-EPS&&o1>=-EPS&&o2>=-EPS){startTri=cur;break;}
                int w=-1;double wv=0;
                if(o0<wv){wv=o0;w=0;} if(o1<wv){wv=o1;w=1;} if(o2<wv){wv=o2;w=2;}
                if(w>=0&&t.neighbors[w]>=0) cur=t.neighbors[w]; else break;
            }
            if(startTri<0){
                for(int i=0;i<(int)tris.size();i++){
                    const auto& t=tris[i];
                    if(t.v[0]<0) continue;
                    if(orient2d(pts[t.v[0]],pts[t.v[1]],p)>=-EPS&&
                       orient2d(pts[t.v[1]],pts[t.v[2]],p)>=-EPS&&
                       orient2d(pts[t.v[2]],pts[t.v[0]],p)>=-EPS){startTri=i;break;}
                }
            }
        }
        if(startTri<0) continue;

        // BFS cavity (epoch-based visited to avoid per-point allocation)
        epoch++;
        visitEpoch.resize(tris.size(), 0);
        std::vector<int> bad;
        std::vector<int> stk={startTri};
        while(!stk.empty()){
            int cur=stk.back(); stk.pop_back();
            if(cur<0||cur>=(int)tris.size()||visitEpoch[cur]==epoch) continue;
            visitEpoch[cur]=epoch;
            const auto& t=tris[cur];
            if(t.v[0]<0) continue;
            if(inCircle(pts[t.v[0]],pts[t.v[1]],pts[t.v[2]],p)>0){
                bad.push_back(cur);
                for(int j=0;j<3;j++){
                    if(t.neighbors[j]>=0 && visitEpoch[t.neighbors[j]]!=epoch)
                        stk.push_back(t.neighbors[j]);
                }
            }
        }

        // Force-insert if cavity empty (near-coincident vertices)
        if(bad.empty()) bad.push_back(startTri);

        // Collect boundary polygon
        std::unordered_set<int> badSet(bad.begin(), bad.end());
        struct BndEdge { int v0, v1, outerTri, outerLocalEdge; };
        std::vector<BndEdge> poly;
        for(int bi : bad){
            const auto& t = tris[bi];
            for(int j=0;j<3;j++){
                int nb=t.neighbors[j];
                if(nb<0||!badSet.count(nb)){
                    int oe=-1;
                    if(nb>=0) for(int k=0;k<3;k++) if(tris[nb].neighbors[k]==bi){oe=k;break;}
                    poly.push_back({t.v[j], t.v[(j+1)%3], nb, oe});
                }
            }
        }

        int nNew=(int)poly.size();
        std::sort(bad.begin(), bad.end());
        std::vector<int> slots;
        for(int i=0; i<nNew && i<(int)bad.size(); i++) slots.push_back(bad[i]);
        while((int)slots.size()<nNew){
            slots.push_back((int)tris.size());
            tris.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
        }
        for(int i=nNew; i<(int)bad.size(); i++){
            tris[bad[i]].v={-1,-1,-1};
            tris[bad[i]].neighbors={-1,-1,-1};
        }

        for(int i=0; i<nNew; i++){
            int s=slots[i];
            tris[s].v = {poly[i].v0, poly[i].v1, pidx};
            tris[s].neighbors = {-1,-1,-1};
            tris[s].region_attrib=0; tris[s].region_max_area=-1; tris[s].marker=0;
            if(orient2d(pts[tris[s].v[0]],pts[tris[s].v[1]],pts[tris[s].v[2]])<0)
                std::swap(tris[s].v[0], tris[s].v[1]);
        }

        // Wire outer adjacency
        for(int i=0; i<nNew; i++){
            int s=slots[i];
            if(poly[i].outerTri<0) continue;
            for(int j=0;j<3;j++){
                int a=tris[s].v[j], b=tris[s].v[(j+1)%3];
                if((a==poly[i].v0&&b==poly[i].v1)||(a==poly[i].v1&&b==poly[i].v0)){
                    tris[s].neighbors[j]=poly[i].outerTri;
                    if(poly[i].outerLocalEdge>=0)
                        tris[poly[i].outerTri].neighbors[poly[i].outerLocalEdge]=s;
                    break;
                }
            }
        }

        // Wire internal adjacency
        std::unordered_map<int, std::pair<int,int>> pidxEdges;
        for(int i=0; i<nNew; i++){
            int s=slots[i];
            for(int j=0;j<3;j++){
                int a=tris[s].v[j], b=tris[s].v[(j+1)%3];
                if(a!=pidx && b!=pidx) continue;
                int other = (a==pidx)?b:a;
                auto it=pidxEdges.find(other);
                if(it!=pidxEdges.end()){
                    auto [os,oe]=it->second;
                    tris[s].neighbors[j]=os;
                    tris[os].neighbors[oe]=s;
                    pidxEdges.erase(it);
                } else {
                    pidxEdges[other]={s,j};
                }
            }
        }
        lastTri = slots[0];
    }

    // Remove dead triangles and super-triangle vertices
    std::vector<bool> keep(tris.size(), false);
    for(int i=0;i<(int)tris.size();i++){
        if(tris[i].v[0]<0) continue;
        if(tris[i].v[0]>=n||tris[i].v[1]>=n||tris[i].v[2]>=n) continue;
        keep[i]=true;
    }
    std::vector<int> remap(tris.size(), -1);
    std::vector<Triangle> newTris;
    for(int i=0;i<(int)tris.size();i++)
        if(keep[i]){remap[i]=(int)newTris.size(); newTris.push_back(tris[i]);}
    for(auto& t:newTris)
        for(int j=0;j<3;j++) t.neighbors[j]=(t.neighbors[j]>=0)?remap[t.neighbors[j]]:-1;
    tris=std::move(newTris);
    pts.resize(n);

    // Detect and re-insert orphaned vertices
    mesh.rebuildAdjacency();
    std::vector<int> orphanVisitEpoch;
    int orphanEpoch = 0;
    for(int retry=0; retry<10; retry++){
        std::vector<bool> inMesh(n, false);
        for(auto& t : tris)
            if(t.v[0]>=0) for(int j=0;j<3;j++) if(t.v[j]<n) inMesh[t.v[j]]=true;

        std::vector<int> orphans;
        for(int i=0;i<n;i++) if(!inMesh[i]) orphans.push_back(i);
        if(orphans.empty()) break;

        for(int pidx : orphans){
            const Point& p = pts[pidx];
            int containTri = mesh.locateTriangle(p.x, p.y);
            if(containTri<0) continue;

            // Skip if nearly coincident with an existing vertex
            bool tooClose = false;
            for(int j=0;j<3;j++){
                double dx2=pts[tris[containTri].v[j]].x-p.x, dy2=pts[tris[containTri].v[j]].y-p.y;
                if(dx2*dx2+dy2*dy2 < EPS*EPS){ tooClose=true; break; }
            }
            if(tooClose) continue;

            // Bowyer-Watson insertion with strict inCircle to avoid over-expanding
            // the cavity and orphaning nearby vertices
            orphanEpoch++;
            orphanVisitEpoch.resize(tris.size(), 0);
            std::vector<int> bad;
            std::vector<int> stk={containTri};
            while(!stk.empty()){
                int cur=stk.back(); stk.pop_back();
                if(cur<0||cur>=(int)tris.size()||orphanVisitEpoch[cur]==orphanEpoch) continue;
                orphanVisitEpoch[cur]=orphanEpoch;
                const auto& t=tris[cur];
                if(t.v[0]<0) continue;
                if(inCircle(pts[t.v[0]],pts[t.v[1]],pts[t.v[2]],p)>0){
                    bad.push_back(cur);
                    for(int j=0;j<3;j++){
                        if(t.neighbors[j]>=0 && orphanVisitEpoch[t.neighbors[j]]!=orphanEpoch)
                            stk.push_back(t.neighbors[j]);
                    }
                }
            }
            if(bad.empty()) bad.push_back(containTri);

            // Collect boundary polygon
            std::unordered_set<int> badSet(bad.begin(), bad.end());
            struct BndEdge { int v0, v1, outerTri, outerLocalEdge; };
            std::vector<BndEdge> poly;
            for(int bi : bad){
                const auto& t = tris[bi];
                for(int j=0;j<3;j++){
                    int nb=t.neighbors[j];
                    if(nb<0||!badSet.count(nb)){
                        int oe=-1;
                        if(nb>=0) for(int k=0;k<3;k++) if(tris[nb].neighbors[k]==bi){oe=k;break;}
                        poly.push_back({t.v[j], t.v[(j+1)%3], nb, oe});
                    }
                }
            }

            int nNew=(int)poly.size();
            std::sort(bad.begin(), bad.end());
            std::vector<int> slots;
            for(int i=0; i<nNew && i<(int)bad.size(); i++) slots.push_back(bad[i]);
            while((int)slots.size()<nNew){
                slots.push_back((int)tris.size());
                tris.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
            }
            for(int i=nNew; i<(int)bad.size(); i++){
                tris[bad[i]].v={-1,-1,-1};
                tris[bad[i]].neighbors={-1,-1,-1};
            }

            for(int i=0; i<nNew; i++){
                int s=slots[i];
                tris[s].v = {poly[i].v0, poly[i].v1, pidx};
                tris[s].neighbors = {-1,-1,-1};
                tris[s].region_attrib=0; tris[s].region_max_area=-1; tris[s].marker=0;
                if(orient2d(pts[tris[s].v[0]],pts[tris[s].v[1]],pts[tris[s].v[2]])<0)
                    std::swap(tris[s].v[0], tris[s].v[1]);
            }

            // Wire outer adjacency
            for(int i=0; i<nNew; i++){
                int s=slots[i];
                if(poly[i].outerTri<0) continue;
                for(int j=0;j<3;j++){
                    int a=tris[s].v[j], b=tris[s].v[(j+1)%3];
                    if((a==poly[i].v0&&b==poly[i].v1)||(a==poly[i].v1&&b==poly[i].v0)){
                        tris[s].neighbors[j]=poly[i].outerTri;
                        if(poly[i].outerLocalEdge>=0)
                            tris[poly[i].outerTri].neighbors[poly[i].outerLocalEdge]=s;
                        break;
                    }
                }
            }

            // Wire internal adjacency
            std::unordered_map<int, std::pair<int,int>> pidxEdges;
            for(int i=0; i<nNew; i++){
                int s=slots[i];
                for(int j=0;j<3;j++){
                    int a=tris[s].v[j], b=tris[s].v[(j+1)%3];
                    if(a!=pidx && b!=pidx) continue;
                    int other = (a==pidx)?b:a;
                    auto it=pidxEdges.find(other);
                    if(it!=pidxEdges.end()){
                        auto [os,oe]=it->second;
                        tris[s].neighbors[j]=os;
                        tris[os].neighbors[oe]=s;
                        pidxEdges.erase(it);
                    } else {
                        pidxEdges[other]={s,j};
                    }
                }
            }
        }
    }
}

// ============================================================
// CDT: enforce PSLG segments via cavity retriangulation
// ============================================================

int enforceConstraints(Mesh& mesh, bool quiet){
    if(mesh.segments.empty()) return 0;

    auto buildV2T = [&]() {
        std::vector<std::vector<int>> v2t(mesh.vertices.size());
        for(int i=0;i<(int)mesh.triangles.size();i++){
            if(mesh.triangles[i].v[0]<0) continue;
            for(int j=0;j<3;j++) v2t[mesh.triangles[i].v[j]].push_back(i);
        }
        return v2t;
    };

    // Ear clipping triangulation of a simple polygon
    auto earClip = [&](const std::vector<int>& poly) -> std::vector<std::array<int,3>> {
        std::vector<std::array<int,3>> result;
        std::vector<int> p = poly;
        while(p.size()>3){
            bool clipped=false;
            int n=(int)p.size();
            // Two passes: first strict, then relaxed tolerance for thin polygons
            for(int pass2=0;pass2<2&&!clipped;pass2++){
                double tol=(pass2==0)?0.0:-1e-8;
                for(int i=0;i<n;i++){
                    int prev=p[(i+n-1)%n], cur=p[i], next=p[(i+1)%n];
                    if(orient2d(mesh.vertices[prev],mesh.vertices[cur],mesh.vertices[next])<=tol) continue;
                    bool inside=false;
                    for(int k=0;k<n&&!inside;k++){
                        if(k==(i+n-1)%n||k==i||k==(i+1)%n) continue;
                        int vi=p[k];
                        // Use strict positive test: vertices on ear boundary are OK
                        if(orient2d(mesh.vertices[prev],mesh.vertices[cur],mesh.vertices[vi])>0&&
                           orient2d(mesh.vertices[cur],mesh.vertices[next],mesh.vertices[vi])>0&&
                           orient2d(mesh.vertices[next],mesh.vertices[prev],mesh.vertices[vi])>0)
                            inside=true;
                    }
                    if(inside) continue;
                    result.push_back({prev,cur,next});
                    p.erase(p.begin()+i);
                    clipped=true;
                    break;
                }
            }
            if(!clipped) break;
        }
        // Handle remaining polygon (3 or more vertices that couldn't be ear-clipped)
        if(p.size()>=3){
            if(p.size()==3){
                double o=orient2d(mesh.vertices[p[0]],mesh.vertices[p[1]],mesh.vertices[p[2]]);
                if(o>0) result.push_back({p[0],p[1],p[2]});
                else if(o<0) result.push_back({p[0],p[2],p[1]});
            } else {
                // Fallback: fan triangulation from first vertex
                for(int i=1;i<(int)p.size()-1;i++){
                    double o=orient2d(mesh.vertices[p[0]],mesh.vertices[p[i]],mesh.vertices[p[i+1]]);
                    if(o>0) result.push_back({p[0],p[i],p[i+1]});
                    else if(o<0) result.push_back({p[0],p[i+1],p[i]});
                    // Skip degenerate (collinear) triangles
                }
            }
        }
        return result;
    };

    auto v2t = buildV2T();

    // Split segments at collinear intermediate vertices (Triangle's approach).
    // When vertices lie exactly on a segment between its endpoints, the walk-based
    // CDT can't find a starting triangle (orient2d returns 0 for all neighbors).
    // Splitting into sub-segments makes each one trivially enforceable.
    {
        // Sort vertex indices by x-coordinate for efficient range queries
        std::vector<int> xorder(mesh.vertices.size());
        std::iota(xorder.begin(), xorder.end(), 0);
        std::sort(xorder.begin(), xorder.end(), [&](int a, int b){
            return mesh.vertices[a].x < mesh.vertices[b].x;
        });
        // Parallel array of sorted x-values for binary search
        std::vector<double> xsorted(xorder.size());
        for(int i=0;i<(int)xorder.size();i++) xsorted[i]=mesh.vertices[xorder[i]].x;

        std::vector<Segment> newSegs;
        for(auto& seg : mesh.segments){
            int su=seg.v0, sv=seg.v1;
            double sx=mesh.vertices[su].x, sy=mesh.vertices[su].y;
            double dx=mesh.vertices[sv].x-sx, dy=mesh.vertices[sv].y-sy;
            double len2=dx*dx+dy*dy;
            if(len2<EPS*EPS){ newSegs.push_back(seg); continue; }

            // Find all vertices strictly between su and sv on the line segment
            double EPS2len2 = EPS * EPS * len2; // squared collinearity threshold
            double bxlo=std::min(sx,sx+dx)-EPS, bxhi=std::max(sx,sx+dx)+EPS;
            double bylo=std::min(sy,sy+dy)-EPS, byhi=std::max(sy,sy+dy)+EPS;
            // Binary search for x-range [bxlo, bxhi]
            int ilo=(int)(std::lower_bound(xsorted.begin(),xsorted.end(),bxlo)-xsorted.begin());
            int ihi=(int)(std::upper_bound(xsorted.begin(),xsorted.end(),bxhi)-xsorted.begin());
            std::vector<std::pair<double,int>> onSeg; // (parameter t, vertex index)
            for(int ii=ilo;ii<ihi;ii++){
                int i=xorder[ii];
                if(i==su||i==sv) continue;
                double vy=mesh.vertices[i].y;
                if(vy<bylo||vy>byhi) continue;
                double px=mesh.vertices[i].x-sx, py=vy-sy;
                // Check collinearity: |cross|² ≤ EPS² * len2
                double cross=dx*py-dy*px;
                if(cross*cross > EPS2len2) continue;
                // Project onto segment to get parameter t
                double t=(dx*px+dy*py)/len2;
                if(t<=EPS_SCALE||t>=1.0-EPS_SCALE) continue;
                onSeg.push_back({t, i});
            }

            if(onSeg.empty()){
                newSegs.push_back(seg);
            } else {
                std::sort(onSeg.begin(), onSeg.end());
                int prev=su;
                for(auto& [t,vi] : onSeg){
                    newSegs.push_back({prev, vi, seg.marker});
                    prev=vi;
                }
                newSegs.push_back({prev, sv, seg.marker});
            }
        }
        if(newSegs.size() != mesh.segments.size()){
            if(!quiet)
                std::cerr<<"  Split "<<mesh.segments.size()<<" segments into "<<newSegs.size()<<" sub-segments (collinear vertices)\n";
            mesh.segments = newSegs;
        }
    }

    int totalMissing = 0;
    for(int pass=0; pass<5; pass++){
        int enforced = 0;
        for(auto& seg : mesh.segments){
            int u=seg.v0, v=seg.v1;
            if(mesh.hasEdge(u,v,v2t)) continue;

            // Try walking from both endpoints, with optional SoS perturbation
            bool success = false;
            for(int attempt=0; attempt<4 && !success; attempt++){
                int su = (attempt%2==0)?u:v;
                int sv = (attempt%2==0)?v:u;

                // For attempts 2-3, use SoS perturbation to break collinearity
                // Create a virtual target point slightly offset from sv
                Point sv_pert = mesh.vertices[sv];
                if(attempt >= 2){
                    // Perpendicular perturbation to the segment direction
                    double dx=mesh.vertices[sv].x-mesh.vertices[su].x;
                    double dy=mesh.vertices[sv].y-mesh.vertices[su].y;
                    double len=std::sqrt(dx*dx+dy*dy);
                    if(len>EPS){
                        // Perturb perpendicular: rotate direction 90° and scale tiny
                        double eps=len*1e-8;
                        sv_pert.y += eps; // simple y perturbation
                    }
                }

                // orient2d test using potentially perturbed sv
                auto orient_seg = [&](const Point& p) -> double {
                    return orient2d(mesh.vertices[su], sv_pert, p);
                };

                // Find starting triangle using finddirection
                int startTri = -1;
                {
                    int curTri = -1;
                    for(int ti : v2t[su])
                        if(mesh.triangles[ti].v[0]>=0){ curTri=ti; break; }
                    if(curTri<0) continue;

                    // Check if sv is already adjacent
                    for(int ti : v2t[su]){
                        auto& t=mesh.triangles[ti];
                        if(t.v[0]<0) continue;
                        for(int j=0;j<3;j++) if(t.v[j]==sv){success=true;break;}
                        if(success) break;
                    }
                    if(success) continue;

                    auto getLocalU = [&](int ti) -> int {
                        for(int j=0;j<3;j++) if(mesh.triangles[ti].v[j]==su) return j;
                        return -1;
                    };

                    // Scan all fan triangles to find the one the segment ray exits through.
                    // For each triangle, check: right is on/above and left is on/below the
                    // line su→sv, AND the opposite edge midpoint is in the forward direction.
                    int bestTri = -1;
                    double bestDot = -1e30;
                    double sdx=mesh.vertices[sv].x-mesh.vertices[su].x;
                    double sdy=mesh.vertices[sv].y-mesh.vertices[su].y;
                    for(int ti : v2t[su]){
                        auto& t=mesh.triangles[ti];
                        if(t.v[0]<0) continue;
                        int lu = getLocalU(ti);
                        if(lu<0) continue;
                        int right = t.v[(lu+1)%3], left = t.v[(lu+2)%3];
                        double oRight = orient_seg(mesh.vertices[right]);
                        double oLeft  = orient_seg(mesh.vertices[left]);
                        if(oRight >= 0 && oLeft <= 0){
                            // Straddling triangle — check forward direction
                            double mx=(mesh.vertices[right].x+mesh.vertices[left].x)/2;
                            double my=(mesh.vertices[right].y+mesh.vertices[left].y)/2;
                            double dmx=mx-mesh.vertices[su].x, dmy=my-mesh.vertices[su].y;
                            double dot=sdx*dmx+sdy*dmy;
                            if(dot > bestDot){
                                bestDot = dot;
                                bestTri = ti;
                            }
                        }
                    }
                    startTri = bestTri;
                }
                if(success) continue;
                if(startTri<0) continue;

                // Walk from su to sv collecting crossed triangles and chains
                auto& st = mesh.triangles[startTri];
                int localU = -1;
                for(int j=0;j<3;j++) if(st.v[j]==su) localU=j;
                int a=st.v[(localU+1)%3], b=st.v[(localU+2)%3];
                double oa=orient_seg(mesh.vertices[a]);

                int firstAbove, firstBelow;
                if(oa>=0){ firstAbove=a; firstBelow=b; }
                else     { firstAbove=b; firstBelow=a; }

                std::vector<int> crossedTris = {startTri};
                std::vector<int> aboveChain = {firstAbove};
                std::vector<int> belowChain = {firstBelow};

                // Collect boundary edges (edges of crossed tris not crossed by segment)
                struct BndEdge { int v0,v1,outerTri,outerEdge; };
                std::vector<BndEdge> bndEdges;

                // Boundary edges of first triangle (the two edges incident to su)
                for(int j=0;j<3;j++){
                    int ev0=st.v[j], ev1=st.v[(j+1)%3];
                    if((ev0==firstAbove&&ev1==firstBelow)||(ev0==firstBelow&&ev1==firstAbove)) continue;
                    int nb=st.neighbors[j]; int oe=-1;
                    if(nb>=0) for(int k=0;k<3;k++) if(mesh.triangles[nb].neighbors[k]==startTri){oe=k;break;}
                    bndEdges.push_back({ev0,ev1,nb,oe});
                }

                // Check if first triangle already contains sv
                bool reachedSv = false;
                for(int j=0;j<3;j++) if(st.v[j]==sv){ reachedSv=true; break; }

                int prevAbove=firstAbove, prevBelow=firstBelow;
                std::set<int> visitedTris;
                visitedTris.insert(startTri);

                while(!reachedSv){
                    int lastTri2 = crossedTris.back();
                    auto& lt = mesh.triangles[lastTri2];
                    // Find the neighbor across the frontier edge (prevAbove, prevBelow)
                    int nextTri = -1;
                    for(int j=0;j<3;j++){
                        int ev0=lt.v[j], ev1=lt.v[(j+1)%3];
                        if((ev0==prevAbove&&ev1==prevBelow)||(ev0==prevBelow&&ev1==prevAbove)){
                            nextTri=lt.neighbors[j]; break;
                        }
                    }
                    if(nextTri<0 || visitedTris.count(nextTri)) break;
                    visitedTris.insert(nextTri);
                    crossedTris.push_back(nextTri);

                    auto& nt = mesh.triangles[nextTri];
                    int newVert = -1;
                    for(int j=0;j<3;j++)
                        if(nt.v[j]!=prevAbove && nt.v[j]!=prevBelow) newVert=nt.v[j];
                    if(newVert<0) break;

                    if(newVert==sv){
                        reachedSv=true;
                        // Add boundary edges of last triangle (not the frontier edge)
                        for(int j=0;j<3;j++){
                            int ev0=nt.v[j], ev1=nt.v[(j+1)%3];
                            if((ev0==prevAbove&&ev1==prevBelow)||(ev0==prevBelow&&ev1==prevAbove)) continue;
                            int nb=nt.neighbors[j]; int oe=-1;
                            if(nb>=0) for(int k=0;k<3;k++) if(mesh.triangles[nb].neighbors[k]==nextTri){oe=k;break;}
                            bndEdges.push_back({ev0,ev1,nb,oe});
                        }
                    } else {
                        double side=orient_seg(mesh.vertices[newVert]);
                        if(side>0){
                            aboveChain.push_back(newVert);
                            // Boundary edge: the edge from newVert to prevAbove
                            for(int j=0;j<3;j++){
                                int ev0=nt.v[j], ev1=nt.v[(j+1)%3];
                                if((ev0==prevAbove&&ev1==prevBelow)||(ev0==prevBelow&&ev1==prevAbove)) continue;
                                if((ev0==newVert&&ev1==prevBelow)||(ev0==prevBelow&&ev1==newVert)) continue;
                                int nb=nt.neighbors[j]; int oe=-1;
                                if(nb>=0) for(int k=0;k<3;k++) if(mesh.triangles[nb].neighbors[k]==nextTri){oe=k;break;}
                                bndEdges.push_back({ev0,ev1,nb,oe});
                            }
                            prevAbove=newVert;
                        } else {
                            belowChain.push_back(newVert);
                            for(int j=0;j<3;j++){
                                int ev0=nt.v[j], ev1=nt.v[(j+1)%3];
                                if((ev0==prevAbove&&ev1==prevBelow)||(ev0==prevBelow&&ev1==prevAbove)) continue;
                                if((ev0==prevAbove&&ev1==newVert)||(ev0==newVert&&ev1==prevAbove)) continue;
                                int nb=nt.neighbors[j]; int oe=-1;
                                if(nb>=0) for(int k=0;k<3;k++) if(mesh.triangles[nb].neighbors[k]==nextTri){oe=k;break;}
                                bndEdges.push_back({ev0,ev1,nb,oe});
                            }
                            prevBelow=newVert;
                        }
                    }
                }

                if(!reachedSv) continue;

                // Build above and below polygons
                std::vector<int> abovePoly = {su};
                abovePoly.insert(abovePoly.end(), aboveChain.begin(), aboveChain.end());
                abovePoly.push_back(sv);

                std::vector<int> belowPoly = {sv};
                for(int i=(int)belowChain.size()-1;i>=0;i--) belowPoly.push_back(belowChain[i]);
                belowPoly.push_back(su);

                auto aboveTris = earClip(abovePoly);
                auto belowTris = earClip(belowPoly);

                int totalNew = (int)aboveTris.size()+(int)belowTris.size();
                int totalOld = (int)crossedTris.size();
                if(totalNew==0) continue;

                // Assign slots
                std::sort(crossedTris.begin(), crossedTris.end());
                std::vector<int> slots;
                for(int i=0;i<totalOld&&i<totalNew;i++) slots.push_back(crossedTris[i]);
                while((int)slots.size()<totalNew){
                    slots.push_back((int)mesh.triangles.size());
                    mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
                }
                for(int i=totalNew;i<totalOld;i++){
                    mesh.triangles[crossedTris[i]].v={-1,-1,-1};
                    mesh.triangles[crossedTris[i]].neighbors={-1,-1,-1};
                }

                // Place new triangles
                std::vector<std::array<int,3>> allNew;
                allNew.insert(allNew.end(), aboveTris.begin(), aboveTris.end());
                allNew.insert(allNew.end(), belowTris.begin(), belowTris.end());
                for(int i=0;i<totalNew;i++){
                    int s=slots[i];
                    mesh.triangles[s].v=allNew[i];
                    mesh.triangles[s].neighbors={-1,-1,-1};
                    mesh.triangles[s].region_attrib=0;
                    mesh.triangles[s].region_max_area=-1;
                    mesh.triangles[s].marker=0;
                }

                // Wire internal adjacency among new triangles
                std::map<std::pair<int,int>,std::pair<int,int>> edgeMap;
                for(int i=0;i<totalNew;i++){
                    int s=slots[i];
                    auto& t=mesh.triangles[s];
                    for(int j=0;j<3;j++){
                        auto key=edgeKey(t.v[j],t.v[(j+1)%3]);
                        auto it=edgeMap.find(key);
                        if(it!=edgeMap.end()){
                            auto [os,oe]=it->second;
                            t.neighbors[j]=os; mesh.triangles[os].neighbors[oe]=s;
                            edgeMap.erase(it);
                        } else edgeMap[key]={s,j};
                    }
                }

                // Wire boundary edges
                for(auto& be : bndEdges){
                    auto key=edgeKey(be.v0,be.v1);
                    auto it=edgeMap.find(key);
                    if(it!=edgeMap.end()){
                        auto [s,j]=it->second;
                        mesh.triangles[s].neighbors[j]=be.outerTri;
                        if(be.outerTri>=0 && be.outerEdge>=0)
                            mesh.triangles[be.outerTri].neighbors[be.outerEdge]=s;
                    }
                }

                v2t = buildV2T();
                success = true;
                enforced++;
            }
        }
        if(enforced==0 && pass>0) break;
    }

    // Direct adjacency check: for short segments where su and sv are in adjacent
    // triangles (separated by a single edge), insert the edge by flipping.
    {
        mesh.rebuildAdjacency();
        v2t = buildV2T();
        bool anyFixed = true;
        while(anyFixed){
            anyFixed = false;
            for(auto& seg : mesh.segments){
                int su2=seg.v0, sv2=seg.v1;
                if(mesh.hasEdge(su2,sv2,v2t)) continue;
                // Check if su2 and sv2 are in neighboring triangles
                for(int ti : v2t[su2]){
                    if(mesh.hasEdge(su2,sv2,v2t)) break;
                    auto& t=mesh.triangles[ti]; if(t.v[0]<0) continue;
                    for(int j=0;j<3;j++){
                        int nb=t.neighbors[j];
                        if(nb<0) continue;
                        auto& tn=mesh.triangles[nb];
                        bool hasSv=false;
                        for(int k=0;k<3;k++) if(tn.v[k]==sv2){hasSv=true;break;}
                        if(!hasSv) continue;
                        // Found: su2 is in ti, sv2 is in nb, they share edge at j
                        // The shared edge is (t.v[j], t.v[(j+1)%3])
                        int a=t.v[j], b=t.v[(j+1)%3];
                        if(a==su2||b==su2||a==sv2||b==sv2) continue; // shared edge includes endpoints
                        // su2 is opposite the shared edge in ti
                        // sv2 is opposite the shared edge in nb
                        // Flip (a,b) -> (su2, sv2)
                        int p=su2, q=sv2;
                        double o1=orient2d(mesh.vertices[p],mesh.vertices[a],mesh.vertices[q]);
                        double o2=orient2d(mesh.vertices[p],mesh.vertices[q],mesh.vertices[b]);
                        if(o1<0||o2<0) continue; // not convex
                        // Do the flip
                        int N_pa=-1, N_bp=-1;
                        for(int k=0;k<3;k++){
                            int e0=t.v[k], e1=t.v[(k+1)%3];
                            if((e0==p&&e1==a)||(e0==a&&e1==p)) N_pa=t.neighbors[k];
                            if((e0==b&&e1==p)||(e0==p&&e1==b)) N_bp=t.neighbors[k];
                        }
                        int N_aq=-1, N_qb=-1;
                        for(int k=0;k<3;k++){
                            int e0=tn.v[k], e1=tn.v[(k+1)%3];
                            if((e0==a&&e1==q)||(e0==q&&e1==a)) N_aq=tn.neighbors[k];
                            if((e0==q&&e1==b)||(e0==b&&e1==q)) N_qb=tn.neighbors[k];
                        }
                        t.v={p,a,q}; t.neighbors={N_pa, N_aq, nb};
                        tn.v={p,q,b}; tn.neighbors={ti, N_qb, N_bp};
                        if(N_aq>=0) for(int k=0;k<3;k++) if(mesh.triangles[N_aq].neighbors[k]==nb){mesh.triangles[N_aq].neighbors[k]=ti;break;}
                        if(N_bp>=0) for(int k=0;k<3;k++) if(mesh.triangles[N_bp].neighbors[k]==ti){mesh.triangles[N_bp].neighbors[k]=nb;break;}
                        v2t=buildV2T();
                        anyFixed=true;
                        goto nextSeg;
                    }
                }
                nextSeg:;
            }
        }
    }
    {
        // Build set of constraint edges to protect during flipping
        EdgeSet segEdges;
        for(auto& s2:mesh.segments)
            if(mesh.hasEdge(s2.v0,s2.v1,v2t)) segEdges.insert(edgeKey(s2.v0,s2.v1));

        // Check if segment (v1,v2) properly intersects edge (v3,v4)
        auto segsCross = [&](int v1, int v2, int v3, int v4) -> bool {
            double d1=orient2d(mesh.vertices[v1],mesh.vertices[v2],mesh.vertices[v3]);
            double d2=orient2d(mesh.vertices[v1],mesh.vertices[v2],mesh.vertices[v4]);
            if(d1*d2>=0) return false;
            double d3=orient2d(mesh.vertices[v3],mesh.vertices[v4],mesh.vertices[v1]);
            double d4=orient2d(mesh.vertices[v3],mesh.vertices[v4],mesh.vertices[v2]);
            if(d3*d4>=0) return false;
            return true;
        };

        for(auto& seg : mesh.segments){
            int su=seg.v0, sv=seg.v1;
            if(mesh.hasEdge(su,sv,v2t)) continue;

            // Walk from su toward sv collecting crossing edges, then flip them
            {
                // Try walking from both endpoints
                for(int attempt=0; attempt<2 && !mesh.hasEdge(su,sv,v2t); attempt++){
                    int wu = (attempt==0)?su:sv;
                    int wv = (attempt==0)?sv:su;

                    // Find starting triangle in wu's fan that straddles the segment
                    double sdx = mesh.vertices[wv].x - mesh.vertices[wu].x;
                    double sdy = mesh.vertices[wv].y - mesh.vertices[wu].y;
                    int startTri = -1;
                    for(int ti : v2t[wu]){
                        auto& t = mesh.triangles[ti];
                        if(t.v[0]<0) continue;
                        int lu = -1;
                        for(int j=0;j<3;j++) if(t.v[j]==wu){ lu=j; break; }
                        if(lu<0) continue;
                        int right = t.v[(lu+1)%3], left = t.v[(lu+2)%3];
                        double oR = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[right]);
                        double oL = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[left]);
                        if(oR >= 0 && oL <= 0){
                            double mx=(mesh.vertices[right].x+mesh.vertices[left].x)*0.5;
                            double my=(mesh.vertices[right].y+mesh.vertices[left].y)*0.5;
                            if(sdx*(mx-mesh.vertices[wu].x)+sdy*(my-mesh.vertices[wu].y) > 0){
                                startTri=ti; break;
                            }
                        }
                    }
                    if(startTri<0) continue;

                    // Walk to collect crossing edges
                    auto& st = mesh.triangles[startTri];
                    int lu = -1;
                    for(int j=0;j<3;j++) if(st.v[j]==wu){ lu=j; break; }
                    int ra=st.v[(lu+1)%3], rb=st.v[(lu+2)%3];
                    double oa = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[ra]);
                    int prevAbove = (oa>=0)?ra:rb, prevBelow = (oa>=0)?rb:ra;

                    std::deque<std::pair<int,int>> flipQueue;
                    flipQueue.push_back(edgeKey(prevAbove, prevBelow));

                    std::unordered_set<int> visited;
                    visited.insert(startTri);
                    int curTri = startTri;
                    bool reached = false;
                    for(int step=0; step<(int)mesh.triangles.size()*2 && !reached; step++){
                        int nextTri = -1;
                        auto& ct = mesh.triangles[curTri];
                        for(int j=0;j<3;j++){
                            int ev0=ct.v[j], ev1=ct.v[(j+1)%3];
                            if((ev0==prevAbove&&ev1==prevBelow)||(ev0==prevBelow&&ev1==prevAbove)){
                                nextTri=ct.neighbors[j]; break;
                            }
                        }
                        if(nextTri<0 || visited.count(nextTri)) break;
                        visited.insert(nextTri);
                        curTri = nextTri;

                        auto& nt = mesh.triangles[nextTri];
                        int newVert = -1;
                        for(int j=0;j<3;j++) if(nt.v[j]!=prevAbove && nt.v[j]!=prevBelow) newVert=nt.v[j];
                        if(newVert<0) break;
                        if(newVert==wv){ reached=true; break; }

                        double side = orient2d(mesh.vertices[wu], mesh.vertices[wv], mesh.vertices[newVert]);
                        if(side > 0){
                            flipQueue.push_back(edgeKey(newVert, prevBelow));
                            prevAbove = newVert;
                        } else {
                            flipQueue.push_back(edgeKey(prevAbove, newVert));
                            prevBelow = newVert;
                        }
                    }
                    if(!reached) continue;

                    // Flip crossing edges (standard CDT flip algorithm)
                    int maxIter = (int)flipQueue.size() * (int)flipQueue.size() + 10;
                    while(!flipQueue.empty() && maxIter-- > 0){
                        auto [ea, eb] = flipQueue.front(); flipQueue.pop_front();
                        if(!segsCross(su, sv, ea, eb)) continue;
                        if(segEdges.count(edgeKey(ea, eb))) continue;

                        // Find two triangles sharing edge (ea,eb) via v2t
                        int ti=-1, nbTri=-1, tiEdge=-1, nbEdge2=-1;
                        for(int t : v2t[ea]){
                            auto& tr = mesh.triangles[t];
                            if(tr.v[0]<0) continue;
                            for(int j=0;j<3;j++){
                                if((tr.v[j]==ea && tr.v[(j+1)%3]==eb)||(tr.v[j]==eb && tr.v[(j+1)%3]==ea)){
                                    if(ti<0){ ti=t; tiEdge=j; }
                                    else { nbTri=t; nbEdge2=j; }
                                }
                            }
                        }
                        if(ti<0 || nbTri<0) continue;

                        int p = mesh.triangles[ti].v[(tiEdge+2)%3];
                        int q = mesh.triangles[nbTri].v[(nbEdge2+2)%3];

                        // Check convexity
                        double co1=orient2d(mesh.vertices[p],mesh.vertices[ea],mesh.vertices[q]);
                        double co2=orient2d(mesh.vertices[p],mesh.vertices[q],mesh.vertices[eb]);
                        if(co1<0||co2<0){
                            flipQueue.push_back({ea, eb}); // defer
                            continue;
                        }

                        // Flip (ea,eb) → (p,q)
                        auto& T1=mesh.triangles[ti];
                        auto& T2=mesh.triangles[nbTri];
                        int N_pa=-1, N_bp=-1;
                        for(int k=0;k<3;k++){
                            int e0=T1.v[k], e1=T1.v[(k+1)%3];
                            if((e0==p&&e1==ea)||(e0==ea&&e1==p)) N_pa=T1.neighbors[k];
                            if((e0==eb&&e1==p)||(e0==p&&e1==eb)) N_bp=T1.neighbors[k];
                        }
                        int N_aq=-1, N_qb=-1;
                        for(int k=0;k<3;k++){
                            int e0=T2.v[k], e1=T2.v[(k+1)%3];
                            if((e0==ea&&e1==q)||(e0==q&&e1==ea)) N_aq=T2.neighbors[k];
                            if((e0==q&&e1==eb)||(e0==eb&&e1==q)) N_qb=T2.neighbors[k];
                        }
                        T1.v={p,ea,q}; T1.neighbors={N_pa, N_aq, nbTri};
                        T2.v={p,q,eb}; T2.neighbors={ti, N_qb, N_bp};
                        if(N_aq>=0) for(int k=0;k<3;k++) if(mesh.triangles[N_aq].neighbors[k]==nbTri){mesh.triangles[N_aq].neighbors[k]=ti;break;}
                        if(N_bp>=0) for(int k=0;k<3;k++) if(mesh.triangles[N_bp].neighbors[k]==ti){mesh.triangles[N_bp].neighbors[k]=nbTri;break;}

                        // Incremental v2t: T1 lost eb gained q, T2 lost ea gained p
                        v2t[eb].erase(std::remove(v2t[eb].begin(), v2t[eb].end(), ti), v2t[eb].end());
                        v2t[ea].erase(std::remove(v2t[ea].begin(), v2t[ea].end(), nbTri), v2t[ea].end());
                        v2t[q].push_back(ti);
                        v2t[p].push_back(nbTri);


                        // If new edge still crosses, re-enqueue
                        if(segsCross(su, sv, p, q))
                            flipQueue.push_back(edgeKey(p, q));
                    }
                }
            }
            if(mesh.hasEdge(su,sv,v2t)){
                segEdges.insert(edgeKey(su,sv));
                continue; // skip cavity approach
            }

            // Second fallback: cavity-based approach
            {
                std::set<int> crossedSet;
                // Include triangles that have at least one edge strictly crossing the segment
                for(int ti=0;ti<(int)mesh.triangles.size();ti++){
                    auto& t=mesh.triangles[ti]; if(t.v[0]<0) continue;
                    for(int j=0;j<3;j++){
                        int a=t.v[j], b=t.v[(j+1)%3];
                        if(segsCross(su,sv,a,b)){
                            crossedSet.insert(ti);
                            break;
                        }
                    }
                }
                // Also include su/sv fan triangles that are "between" the crossed region
                // A fan triangle of su is inside if its opposite edge has one vertex above
                // and one below the segment line
                for(int vi : {su, sv}){
                    for(int ti : v2t[vi]){
                        auto& t=mesh.triangles[ti]; if(t.v[0]<0) continue;
                        // Check if any neighbor of this triangle is in crossedSet
                        for(int j=0;j<3;j++){
                            if(crossedSet.count(t.neighbors[j])){
                                // Check: is there a crossing edge in this triangle?
                                for(int k=0;k<3;k++){
                                    int a=t.v[k], b=t.v[(k+1)%3];
                                    if(segsCross(su,sv,a,b)){crossedSet.insert(ti); break;}
                                }
                                break;
                            }
                        }
                    }
                }

                if(!crossedSet.empty()){
                    // Extract boundary edges of crossed region
                    struct BndEdge { int v0,v1,outerTri; };
                    std::vector<BndEdge> bndEdges;
                    for(int ti : crossedSet){
                        auto& t=mesh.triangles[ti];
                        for(int j=0;j<3;j++){
                            int nb=t.neighbors[j];
                            if(nb<0 || !crossedSet.count(nb)){
                                bndEdges.push_back({t.v[j],t.v[(j+1)%3],nb});
                            }
                        }
                    }

                    // Chain boundary edges into a closed polygon using angular ordering
                    // at vertices with multiple outgoing edges (pinch points)
                    std::map<int,std::vector<std::pair<int,int>>> adjMap; // v -> [(next_v, edge_idx)]
                    for(int i=0;i<(int)bndEdges.size();i++)
                        adjMap[bndEdges[i].v0].push_back({bndEdges[i].v1, i});

                    // At each vertex, sort outgoing edges by angle for consistent traversal
                    for(auto& [v, edges] : adjMap){
                        if(edges.size()>1){
                            std::sort(edges.begin(), edges.end(), [&](auto& a, auto& b){
                                double ax=mesh.vertices[a.first].x-mesh.vertices[v].x;
                                double ay=mesh.vertices[a.first].y-mesh.vertices[v].y;
                                double bx=mesh.vertices[b.first].x-mesh.vertices[v].x;
                                double by=mesh.vertices[b.first].y-mesh.vertices[v].y;
                                return std::atan2(ay,ax) < std::atan2(by,bx);
                            });
                        }
                    }

                    // Chain: at each vertex with multiple outgoing edges, pick the one
                    // that turns most to the right (smallest CW angle from incoming direction)
                    std::vector<int> poly;
                    std::set<int> usedEdges;
                    int startV = bndEdges[0].v0;
                    int curV = startV;
                    int prevV = -1;
                    bool chainOK = true;

                    for(int iter=0;iter<(int)bndEdges.size()+1;iter++){
                        poly.push_back(curV);
                        auto& outs = adjMap[curV];
                        int bestIdx = -1;

                        if(outs.size()==1 && !usedEdges.count(outs[0].second)){
                            bestIdx = 0;
                        } else if(prevV>=0) {
                            // Pick the outgoing edge that makes the smallest CW turn
                            // from the incoming direction (prevV -> curV)
                            double inDx=mesh.vertices[curV].x-mesh.vertices[prevV].x;
                            double inDy=mesh.vertices[curV].y-mesh.vertices[prevV].y;
                            double inAngle=std::atan2(inDy,inDx);
                            double bestTurn=1e30;
                            for(int k=0;k<(int)outs.size();k++){
                                if(usedEdges.count(outs[k].second)) continue;
                                double ox=mesh.vertices[outs[k].first].x-mesh.vertices[curV].x;
                                double oy=mesh.vertices[outs[k].first].y-mesh.vertices[curV].y;
                                double outAngle=std::atan2(oy,ox);
                                // CW turn = how much we turn right (negative = right)
                                double turn=outAngle-inAngle;
                                // Normalize: we want the most-right (most negative) turn
                                // but wraparound: turn in (-pi, pi] mapped to (0, 2pi]
                                // Actually we want to walk the boundary CCW, so pick
                                // leftmost (most CCW) unused edge
                                while(turn<=0) turn+=2*M_PI;
                                if(turn<bestTurn){ bestTurn=turn; bestIdx=k; }
                            }
                        } else {
                            // First step: pick first unused edge
                            for(int k=0;k<(int)outs.size();k++){
                                if(!usedEdges.count(outs[k].second)){ bestIdx=k; break; }
                            }
                        }

                        if(bestIdx<0){ chainOK=false; break; }
                        usedEdges.insert(outs[bestIdx].second);
                        prevV=curV;
                        curV=outs[bestIdx].first;
                        if(curV==startV && (int)poly.size()>2) break;
                    }

                    if(poly.size()>=3 && curV!=startV) chainOK=false;

                    std::vector<std::array<int,3>> newTris;
                    if(chainOK && poly.size()>=3){
                        // Find su and sv in the polygon
                        int suIdx=-1, svIdx=-1;
                        for(int i=0;i<(int)poly.size();i++){
                            if(poly[i]==su) suIdx=i;
                            if(poly[i]==sv) svIdx=i;
                        }

                        if(suIdx>=0 && svIdx>=0){
                            // Split into two polygons at su-sv
                            int n=(int)poly.size();
                            std::vector<int> poly1, poly2;
                            for(int i=suIdx;;){
                                poly1.push_back(poly[i]);
                                if(i==svIdx) break;
                                i=(i+1)%n;
                            }
                            for(int i=svIdx;;){
                                poly2.push_back(poly[i]);
                                if(i==suIdx) break;
                                i=(i+1)%n;
                            }
                            auto t1=earClip(poly1);
                            auto t2=earClip(poly2);
                            newTris.insert(newTris.end(),t1.begin(),t1.end());
                            newTris.insert(newTris.end(),t2.begin(),t2.end());
                        } else {
                            newTris=earClip(poly);
                        }
                    }

                    if(!newTris.empty()){

                        // Assign slots
                        std::vector<int> crossedVec(crossedSet.begin(),crossedSet.end());
                        std::sort(crossedVec.begin(),crossedVec.end());
                        int totalNew=(int)newTris.size();
                        int totalOld=(int)crossedVec.size();

                        std::vector<int> slots;
                        for(int i=0;i<totalOld&&i<totalNew;i++) slots.push_back(crossedVec[i]);
                        while((int)slots.size()<totalNew){
                            slots.push_back((int)mesh.triangles.size());
                            mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
                        }
                        for(int i=totalNew;i<totalOld;i++){
                            mesh.triangles[crossedVec[i]].v={-1,-1,-1};
                            mesh.triangles[crossedVec[i]].neighbors={-1,-1,-1};
                        }

                        for(int i=0;i<totalNew;i++){
                            int s=slots[i];
                            mesh.triangles[s].v=newTris[i];
                            mesh.triangles[s].neighbors={-1,-1,-1};
                        }

                        // Wire internal adjacency
                        std::map<std::pair<int,int>,std::pair<int,int>> edgeMap2;
                        for(int i=0;i<totalNew;i++){
                            int s=slots[i];
                            auto& t=mesh.triangles[s];
                            for(int j=0;j<3;j++){
                                auto key=edgeKey(t.v[j],t.v[(j+1)%3]);
                                auto it=edgeMap2.find(key);
                                if(it!=edgeMap2.end()){
                                    auto [os,oe]=it->second;
                                    t.neighbors[j]=os; mesh.triangles[os].neighbors[oe]=s;
                                    edgeMap2.erase(it);
                                } else edgeMap2[key]={s,j};
                            }
                        }

                        // Wire boundary edges
                        for(auto& be : bndEdges){
                            auto key=edgeKey(be.v0,be.v1);
                            auto it=edgeMap2.find(key);
                            if(it!=edgeMap2.end()){
                                auto [s,j]=it->second;
                                mesh.triangles[s].neighbors[j]=be.outerTri;
                                if(be.outerTri>=0){
                                    auto& ot=mesh.triangles[be.outerTri];
                                    for(int k=0;k<3;k++){
                                        if(crossedSet.count(ot.neighbors[k])){
                                            int e0=ot.v[k], e1=ot.v[(k+1)%3];
                                            if(edgeKey(e0,e1)==key){
                                                ot.neighbors[k]=s;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        v2t=buildV2T();
                    }
                }
            }
            if(mesh.hasEdge(su,sv,v2t)){
                segEdges.insert(edgeKey(su,sv));
            }
        }
    }

    // Count and return remaining missing segments
    int miss=0;
    v2t = buildV2T();
    for(auto& s:mesh.segments){
        if(!mesh.hasEdge(s.v0,s.v1,v2t)){
            miss++;
            if(!quiet)
                std::cerr<<"  UNENFORCED: v"<<s.v0<<" ("<<mesh.vertices[s.v0].x<<","<<mesh.vertices[s.v0].y
                         <<")->v"<<s.v1<<" ("<<mesh.vertices[s.v1].x<<","<<mesh.vertices[s.v1].y
                         <<") len="<<std::sqrt((mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)*(mesh.vertices[s.v1].x-mesh.vertices[s.v0].x)+
                                              (mesh.vertices[s.v1].y-mesh.vertices[s.v0].y)*(mesh.vertices[s.v1].y-mesh.vertices[s.v0].y))<<"\n";
        }
    }
    if(!quiet && miss>0)
        std::cerr<<"Warning: "<<miss<<" of "<<mesh.segments.size()<<" segments could not be enforced.\n";
    return miss;
}

// ============================================================
// Flood-fill hole and exterior removal (bounded by segment edges)
// ============================================================

void removeHoles(Mesh& mesh, const Options& opts){
    if(opts.no_holes) return;
    mesh.rebuildAdjacency();

    EdgeSet segEdges;
    for(auto& s : mesh.segments) segEdges.insert(edgeKey(s.v0,s.v1));

    std::vector<bool> remove(mesh.triangles.size(), false);

    // Remove exterior: flood-fill from boundary triangles whose boundary
    // edge is not a constrained segment.  This removes all triangles
    // outside the PSLG, even when no explicit holes are defined.
    {
        std::queue<int> q;
        for(int i=0;i<(int)mesh.triangles.size();i++){
            if(remove[i]) continue;
            auto& t = mesh.triangles[i];
            for(int j=0;j<3;j++){
                if(t.neighbors[j]>=0) continue;              // not a boundary edge
                auto key = edgeKey(t.v[j], t.v[(j+1)%3]);
                if(segEdges.count(key)) continue;             // boundary is a segment – interior
                if(!remove[i]){ remove[i]=true; q.push(i); }
                break;
            }
        }
        while(!q.empty()){
            int ti = q.front(); q.pop();
            auto& t = mesh.triangles[ti];
            for(int j=0;j<3;j++){
                int nb = t.neighbors[j];
                if(nb<0 || remove[nb]) continue;
                auto key = edgeKey(t.v[j], t.v[(j+1)%3]);
                if(segEdges.count(key)) continue;
                remove[nb] = true;
                q.push(nb);
            }
        }
    }

    // Remove explicit holes via flood fill from seed points
    for(auto& hole : mesh.holes){
        int seedTri = mesh.locateTriangle(hole.x, hole.y);
        if(seedTri<0) continue;

        std::queue<int> q;
        q.push(seedTri);
        remove[seedTri] = true;
        while(!q.empty()){
            int ti = q.front(); q.pop();
            auto& t = mesh.triangles[ti];
            for(int j=0;j<3;j++){
                int nb = t.neighbors[j];
                if(nb<0 || remove[nb]) continue;
                auto key = edgeKey(t.v[j], t.v[(j+1)%3]);
                if(segEdges.count(key)) continue;
                remove[nb] = true;
                q.push(nb);
            }
        }
    }

    std::vector<Triangle> kept;
    std::vector<int> remap(mesh.triangles.size(), -1);
    for(int i=0;i<(int)mesh.triangles.size();i++)
        if(!remove[i]){remap[i]=(int)kept.size(); kept.push_back(mesh.triangles[i]);}
    for(auto& t:kept)
        for(int j=0;j<3;j++) t.neighbors[j]=(t.neighbors[j]>=0)?remap[t.neighbors[j]]:-1;
    mesh.triangles=kept;
}

// ============================================================
// Flood-fill region assignment (bounded by segment edges)
// ============================================================

void assignRegions(Mesh& mesh, const Options& opts){
    if(!opts.regions || mesh.regions.empty()) return;
    mesh.rebuildAdjacency();

    EdgeSet segEdges;
    for(auto& s : mesh.segments) segEdges.insert(edgeKey(s.v0,s.v1));

    for(auto& t : mesh.triangles){ t.region_attrib=0; t.region_max_area=-1; }

    for(auto& reg : mesh.regions){
        int seedTri = mesh.locateTriangle(reg.x, reg.y);
        if(seedTri<0) continue;

        std::queue<int> q;
        std::vector<bool> visited(mesh.triangles.size(), false);
        q.push(seedTri);
        visited[seedTri] = true;
        while(!q.empty()){
            int ti = q.front(); q.pop();
            auto& t = mesh.triangles[ti];
            t.region_attrib = reg.attrib;
            t.region_max_area = reg.max_area;
            for(int j=0;j<3;j++){
                int nb = t.neighbors[j];
                if(nb<0 || visited[nb]) continue;
                auto key = edgeKey(t.v[j], t.v[(j+1)%3]);
                if(segEdges.count(key)) continue;
                visited[nb] = true;
                q.push(nb);
            }
        }
    }
}

// ============================================================
// Point insertion for quality refinement (Lawson flip-based)
// ============================================================

// Insert point into mesh using Lawson flips. Returns new vertex index, or -1 on failure.
// Tracks all affected triangle indices in 'affected' for later re-checking.
// Optional flipBuf avoids per-call heap allocation for the flip stack.
int insertPointLawson(Mesh& mesh, const Point& np, int hintTri,
                      const EdgeSet& constrainedEdges,
                      std::vector<int>& affected,
                      std::vector<std::pair<int,int>>* flipBuf=nullptr){
    int pidx=(int)mesh.vertices.size();
    mesh.vertices.push_back(np);

    int ti=mesh.locateTriangle(np.x, np.y, hintTri);
    if(ti<0){ mesh.vertices.pop_back(); return -1; }

    auto& t=mesh.triangles[ti];
    // Check near-coincident with existing vertex
    for(int j=0;j<3;j++){
        double dx=mesh.vertices[t.v[j]].x-np.x, dy=mesh.vertices[t.v[j]].y-np.y;
        if(dx*dx+dy*dy < EPS*EPS){ mesh.vertices.pop_back(); return -1; }
    }

    // Inherit region from the containing triangle
    double reg_attr=t.region_attrib;
    double reg_area=t.region_max_area;

    // Check if point is on an edge of the containing triangle
    int onEdge=-1;
    for(int j=0;j<3;j++){
        int va=t.v[j], vb=t.v[(j+1)%3];
        double edx=mesh.vertices[vb].x-mesh.vertices[va].x, edy=mesh.vertices[vb].y-mesh.vertices[va].y;
        double o=orient2d_xy(mesh.vertices[va].x,mesh.vertices[va].y,
                             mesh.vertices[vb].x,mesh.vertices[vb].y,np.x,np.y);
        if(o*o < EPS*EPS*(edx*edx+edy*edy)){ onEdge=j; break; }
    }

    std::vector<std::pair<int,int>> flipStackLocal;
    auto& flipStack = flipBuf ? *flipBuf : flipStackLocal;
    flipStack.clear();

    if(onEdge<0){
        // Interior: split triangle into 3
        int v0=t.v[0], v1=t.v[1], v2=t.v[2];
        int n0=t.neighbors[0], n1=t.neighbors[1], n2=t.neighbors[2];

        int t0=ti;
        int t1=(int)mesh.triangles.size();
        mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
        int t2=(int)mesh.triangles.size();
        mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});

        // t0=(v0,v1,pidx), t1=(v1,v2,pidx), t2=(v2,v0,pidx)
        mesh.triangles[t0]={{v0,v1,pidx},{n0,t1,t2},reg_attr,reg_area,0};
        mesh.triangles[t1]={{v1,v2,pidx},{n1,t2,t0},reg_attr,reg_area,0};
        mesh.triangles[t2]={{v2,v0,pidx},{n2,t0,t1},reg_attr,reg_area,0};

        // Fix external back-references (n0 still points to t0=ti, OK)
        if(n1>=0) for(int k=0;k<3;k++)
            if(mesh.triangles[n1].neighbors[k]==ti){ mesh.triangles[n1].neighbors[k]=t1; break; }
        if(n2>=0) for(int k=0;k<3;k++)
            if(mesh.triangles[n2].neighbors[k]==ti){ mesh.triangles[n2].neighbors[k]=t2; break; }

        flipStack.push_back({t0,0});
        flipStack.push_back({t1,0});
        flipStack.push_back({t2,0});
        affected.push_back(t0);
        affected.push_back(t1);
        affected.push_back(t2);
    } else {
        // On edge: split 2 triangles into 4
        int va=t.v[onEdge], vb=t.v[(onEdge+1)%3], vc=t.v[(onEdge+2)%3];
        int n_ab=t.neighbors[onEdge];      // across edge va->vb (where point lies)
        int n_bvc=t.neighbors[(onEdge+1)%3];
        int n_cva=t.neighbors[(onEdge+2)%3];

        if(n_ab<0){
            // Boundary edge: split 1 triangle into 2
            int t0=ti;
            int t1=(int)mesh.triangles.size();
            mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});

            mesh.triangles[t0]={{va,pidx,vc},{-1,t1,n_cva},reg_attr,reg_area,0};
            mesh.triangles[t1]={{pidx,vb,vc},{-1,n_bvc,t0},reg_attr,reg_area,0};

            if(n_bvc>=0) for(int k=0;k<3;k++)
                if(mesh.triangles[n_bvc].neighbors[k]==ti){ mesh.triangles[n_bvc].neighbors[k]=t1; break; }

            flipStack.push_back({t0,2}); // edge vc->va opposite pidx
            flipStack.push_back({t1,1}); // edge vb->vc opposite pidx
            affected.push_back(t0);
            affected.push_back(t1);
        } else {
            // Both sides: split 2 triangles into 4
            auto& tn=mesh.triangles[n_ab];
            int je=-1;
            for(int k=0;k<3;k++) if(tn.neighbors[k]==ti){ je=k; break; }
            if(je<0){ mesh.vertices.pop_back(); return -1; }

            int vd=tn.v[(je+2)%3];
            int n_avd, n_dvb;
            // Determine orientation of shared edge in tn
            if(tn.v[je]==vb && tn.v[(je+1)%3]==va){
                n_avd=tn.neighbors[(je+1)%3]; // across va->vd
                n_dvb=tn.neighbors[(je+2)%3]; // across vd->vb
            } else {
                n_dvb=tn.neighbors[(je+1)%3]; // across vb->vd... = b->q
                n_avd=tn.neighbors[(je+2)%3]; // across vd->va... = q->a
            }

            // Inherit region from the neighbor triangle
            double reg_attr2=tn.region_attrib;
            double reg_area2=tn.region_max_area;

            int t0=ti;
            int t1=(int)mesh.triangles.size();
            mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});
            int t2=n_ab;
            int t3=(int)mesh.triangles.size();
            mesh.triangles.push_back(Triangle{{-1,-1,-1},{-1,-1,-1},0,-1,0});

            mesh.triangles[t0]={{va,pidx,vc},{t3,t1,n_cva},reg_attr,reg_area,0};
            mesh.triangles[t1]={{pidx,vb,vc},{t2,n_bvc,t0},reg_attr,reg_area,0};
            mesh.triangles[t2]={{vb,pidx,vd},{t1,t3,n_dvb},reg_attr2,reg_area2,0};
            mesh.triangles[t3]={{pidx,va,vd},{t0,n_avd,t2},reg_attr2,reg_area2,0};

            // Fix external back-references
            if(n_bvc>=0) for(int k=0;k<3;k++)
                if(mesh.triangles[n_bvc].neighbors[k]==ti){ mesh.triangles[n_bvc].neighbors[k]=t1; break; }
            if(n_avd>=0) for(int k=0;k<3;k++)
                if(mesh.triangles[n_avd].neighbors[k]==n_ab){ mesh.triangles[n_avd].neighbors[k]=t3; break; }

            flipStack.push_back({t0,2}); // vc->va opposite pidx
            flipStack.push_back({t1,1}); // vb->vc opposite pidx
            flipStack.push_back({t2,2}); // vd->vb opposite pidx
            flipStack.push_back({t3,1}); // va->vd opposite pidx
            affected.push_back(t0);
            affected.push_back(t1);
            affected.push_back(t2);
            affected.push_back(t3);
        }
    }

    // Lawson flip loop: restore Delaunay property locally
    while(!flipStack.empty()){
        auto [ft,fj]=flipStack.back();
        flipStack.pop_back();
        auto& tri=mesh.triangles[ft];
        if(tri.v[0]<0) continue;
        int fn=tri.neighbors[fj];
        if(fn<0) continue;

        int fa=tri.v[fj], fb=tri.v[(fj+1)%3], fp=tri.v[(fj+2)%3];

        // Don't flip constrained edges
        if(constrainedEdges.count(edgeKey(fa,fb))) continue;

        auto& tri_n=mesh.triangles[fn];
        if(tri_n.v[0]<0) continue;
        int fj2=-1;
        for(int k=0;k<3;k++) if(tri_n.neighbors[k]==ft){ fj2=k; break; }
        if(fj2<0) continue;

        int fq=tri_n.v[(fj2+2)%3];

        // inCircle test: is fq inside circumcircle of (fa,fb,fp)?
        if(inCircle(mesh.vertices[fa],mesh.vertices[fb],mesh.vertices[fp],mesh.vertices[fq])<=0)
            continue;

        // Check convexity: both new triangles must be CCW
        if(orient2d(mesh.vertices[fa],mesh.vertices[fq],mesh.vertices[fp])<=0) continue;
        if(orient2d(mesh.vertices[fb],mesh.vertices[fp],mesh.vertices[fq])<=0) continue;

        // Save external neighbors before modifying
        int n_fb_fp=tri.neighbors[(fj+1)%3];
        int n_fp_fa=tri.neighbors[(fj+2)%3];
        int n_ext1=tri_n.neighbors[(fj2+1)%3]; // across fa->fq (or fq side)
        int n_ext2=tri_n.neighbors[(fj2+2)%3]; // across fq->fb (or fb side)

        // Verify shared edge orientation and adjust n_ext1/n_ext2 if needed
        if(tri_n.v[fj2]!=fb){
            // tri_n.v[fj2]=fa, tri_n.v[(fj2+1)%3]=fb — reversed from expected
            std::swap(n_ext1, n_ext2);
        }

        double ra=tri.region_attrib, ra2=tri.region_max_area;
        double rb=tri_n.region_attrib, rb2=tri_n.region_max_area;

        // Perform flip: ft=(fa,fq,fp), fn=(fb,fp,fq)
        mesh.triangles[ft]={{fa,fq,fp},{n_ext1,fn,n_fp_fa},ra,ra2,0};
        mesh.triangles[fn]={{fb,fp,fq},{n_fb_fp,ft,n_ext2},rb,rb2,0};

        // Fix external back-references
        // n_ext1 was pointing to fn, now should point to ft
        if(n_ext1>=0) for(int k=0;k<3;k++)
            if(mesh.triangles[n_ext1].neighbors[k]==fn){ mesh.triangles[n_ext1].neighbors[k]=ft; break; }
        // n_fb_fp was pointing to ft, now should point to fn
        if(n_fb_fp>=0) for(int k=0;k<3;k++)
            if(mesh.triangles[n_fb_fp].neighbors[k]==ft){ mesh.triangles[n_fb_fp].neighbors[k]=fn; break; }

        affected.push_back(ft);
        affected.push_back(fn);

        // Push new edges opposite fp for further checking
        // In ft=(fa,fq,fp): fp at idx 2, opposite edge 0: fa->fq
        flipStack.push_back({ft,0});
        // In fn=(fb,fp,fq): fp at idx 1, opposite edge 2: fq->fb
        flipStack.push_back({fn,2});
    }

    // Bump generation for all affected triangles so stale PQ entries are skipped
    for(int ai : affected)
        if(ai>=0 && ai<(int)mesh.triangles.size())
            mesh.triangles[ai].generation++;

    return pidx;
}

// ============================================================
// Quality refinement
// ============================================================

struct TriMetricResult {
    double metric;
    bool angViol;
    bool areaViol;
    double area;
    double effectiveMaxArea;
    double la2, lb2, lc2; // squared edge lengths opposite vertices a, b, c
    double cross;          // 2*signed area = (b-a)x(c-a)
};

TriMetricResult triMetric(const Mesh& mesh, int ti, double minAng, double globalMaxA, bool useRegionArea,
                 double cosMinAng = -1.0){
    auto& t = mesh.triangles[ti];
    if(t.v[0]<0) return {0,false,false,0,-1,0,0,0,0};
    auto& a=mesh.vertices[t.v[0]];
    auto& b=mesh.vertices[t.v[1]];
    auto& c=mesh.vertices[t.v[2]];

    // Squared edge lengths (computed once, reused for both area and angle checks)
    double la2=(b.x-c.x)*(b.x-c.x)+(b.y-c.y)*(b.y-c.y);
    double lb2=(a.x-c.x)*(a.x-c.x)+(a.y-c.y)*(a.y-c.y);
    double lc2=(a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);

    // Area via cross product (no orient2d exact-arithmetic overhead)
    double cross=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    double area=0.5*std::abs(cross);

    double effectiveMaxArea=-1;
    if(useRegionArea && t.region_max_area>0) effectiveMaxArea=t.region_max_area;
    if(globalMaxA>0){
        if(effectiveMaxArea>0) effectiveMaxArea=std::min(effectiveMaxArea,globalMaxA);
        else effectiveMaxArea=globalMaxA;
    }
    // LFS constraint: if any vertex has lfs defined, impose area <= sqrt(3)/4 * lfs^2
    for(int j=0;j<3;j++){
        double vl=mesh.vertices[t.v[j]].lfs;
        if(vl>0){
            constexpr double sqrt3_4=0.43301270189221932; // sqrt(3)/4
            double lfsArea=sqrt3_4*vl*vl;
            if(effectiveMaxArea>0) effectiveMaxArea=std::min(effectiveMaxArea,lfsArea);
            else effectiveMaxArea=lfsArea;
        }
    }
    bool areaViol=(effectiveMaxArea>0 && area>effectiveMaxArea+EPS);

    // Fast angle check using cosine comparison (no hypot/acos)
    // Minimum angle < threshold iff maximum cosine > cosMinAng.
    // cos(angle_opposite_la) = (lb²+lc²-la²) / (2·sqrt(lb²·lc²))
    // Compare squared: (lb²+lc²-la²)² > cos²Threshold · 4·lb²·lc²
    bool angViol=false;
    double maxCosRatio=0; // tracks worst angle violation for metric
    if(minAng>0 && cosMinAng>-0.5){
        double cos2t=cosMinAng*cosMinAng;
        double num, ratio;
        num=lb2+lc2-la2; if(num>0){ ratio=num*num/(cos2t*4.0*lb2*lc2); if(ratio>1.0){ angViol=true; if(ratio>maxCosRatio) maxCosRatio=ratio; }}
        num=la2+lc2-lb2; if(num>0){ ratio=num*num/(cos2t*4.0*la2*lc2); if(ratio>1.0){ angViol=true; if(ratio>maxCosRatio) maxCosRatio=ratio; }}
        num=la2+lb2-lc2; if(num>0){ ratio=num*num/(cos2t*4.0*la2*lb2); if(ratio>1.0){ angViol=true; if(ratio>maxCosRatio) maxCosRatio=ratio; }}
    } else if(minAng>0){
        double ang=minAngle(a,b,c);
        angViol=(ang<minAng-EPS);
        if(angViol) maxCosRatio=(minAng-ang)/minAng+1.0;
    }

    if(!areaViol&&!angViol) return {0,false,false,area,effectiveMaxArea,la2,lb2,lc2,cross};
    double metric=0;
    if(areaViol) metric=area/(effectiveMaxArea>0?effectiveMaxArea:1.0);
    if(angViol) metric=std::max(metric, maxCosRatio);
    return {metric, angViol, areaViol, area, effectiveMaxArea,la2,lb2,lc2,cross};
}

void refineQuality(Mesh& mesh, const Options& opts, int nInputVerts){
    if(!opts.quality&&!opts.area_limit) return;
    // Overshoot the target angle slightly to account for floating-point precision
    // in off-center placement, ensuring we meet the user's requested angle.
    double minAng=opts.quality?opts.min_angle+0.1:0.0;
    double cosMinAng=(minAng>0)?std::cos(minAng*M_PI/180.0):-1.0;
    double globalMaxA=opts.area_limit?opts.max_area:-1.0;
    bool useRegionArea=opts.area_limit;


    // Estimate final triangle count from CDT by computing a mesh size at each
    // vertex, then integrating over each CDT triangle.  Mesh size at a vertex is
    // the minimum of:
    //   (a) area constraint: h = (3/2)*sqrt(amax) where amax is the smallest
    //       max-area from any triangle touching this vertex
    //   (b) local feature size (LFS): distance to nearest non-adjacent segment
    //       or non-neighboring vertex
    // The estimated number of equilateral triangles of side h that fit in an
    // area A is  A / ((sqrt(3)/4)*h^2) = 4*A / (sqrt(3)*h^2).

    int nv=(int)mesh.vertices.size();
    int nt=(int)mesh.triangles.size();

    // Build vertex-to-triangle adjacency
    std::vector<std::vector<int>> v2t(nv);
    for(int i=0;i<nt;i++){
        auto& t=mesh.triangles[i];
        if(t.v[0]<0) continue;
        for(int j=0;j<3;j++) v2t[t.v[j]].push_back(i);
    }

    // Bounding box diagonal — used to cap degenerate circumradii.
    double bbDiag=0;
    {
        double xmin=1e30,xmax=-1e30,ymin=1e30,ymax=-1e30;
        for(int vi=0;vi<nv;vi++){
            xmin=std::min(xmin,mesh.vertices[vi].x);
            xmax=std::max(xmax,mesh.vertices[vi].x);
            ymin=std::min(ymin,mesh.vertices[vi].y);
            ymax=std::max(ymax,mesh.vertices[vi].y);
        }
        bbDiag=std::sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin));
    }

    // Compute local feature size (LFS) at each vertex as min circumradius
    // of incident CDT triangles. Cap at bbox diagonal to ignore degenerate
    // near-zero-area triangles with astronomical circumradii.
    // Also incorporate user-assigned LFS constraints (Point::lfs).
    std::vector<double> lfs(nv, 1e30);
    for(int i=0;i<nt;i++){
        auto& t=mesh.triangles[i];
        if(t.v[0]<0) continue;
        auto& p0=mesh.vertices[t.v[0]];
        auto& p1=mesh.vertices[t.v[1]];
        auto& p2=mesh.vertices[t.v[2]];
        double area=triArea(p0,p1,p2);
        double e0=std::sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
        double e1=std::sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y));
        double e2=std::sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y));
        double R=(area>0)? e0*e1*e2/(4.0*area) : 0;
        if(R>0){
            R=std::min(R, bbDiag);
            for(int j=0;j<3;j++)
                lfs[t.v[j]]=std::min(lfs[t.v[j]], R);
        }
    }
    // Fold in user-assigned LFS constraints
    for(int vi=0;vi<nv;vi++){
        if(mesh.vertices[vi].lfs>0)
            lfs[vi]=std::min(lfs[vi], mesh.vertices[vi].lfs);
    }

    // Compute background mesh size H at each vertex.
    // H = the "local scale" without refinement seeds, from:
    //   (a) Area constraints, and
    //   (b) Median circumradius of incident CDT triangles (captures local geometry).
    // A vertex is a "refinement seed" if its lfs (min circumR) << H (median circumR).
    std::vector<double> bgSize(nv, 1e30);
    for(int vi=0;vi<nv;vi++){
        // Area constraint contribution (only if area limiting is active)
        double amax=globalMaxA;
        if(useRegionArea){
            for(int ti:v2t[vi]){
                auto& t=mesh.triangles[ti];
                if(t.v[0]<0) continue;
                if(t.region_max_area>0)
                    amax=(amax>0)?std::min(amax,t.region_max_area):t.region_max_area;
            }
        }
        if(amax>0) bgSize[vi]=std::sqrt(4.0*amax/std::sqrt(3.0));

        // Median circumradius of incident CDT triangles as local scale
        std::vector<double> incR;
        for(int ti:v2t[vi]){
            auto& t=mesh.triangles[ti];
            if(t.v[0]<0) continue;
            auto& p0=mesh.vertices[t.v[0]];
            auto& p1=mesh.vertices[t.v[1]];
            auto& p2=mesh.vertices[t.v[2]];
            double a2=triArea(p0,p1,p2);
            double ee0=std::sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
            double ee1=std::sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y));
            double ee2=std::sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y));
            double R2=(a2>0)? ee0*ee1*ee2/(4.0*a2) : 0;
            if(R2>0) incR.push_back(std::min(R2, bbDiag));
        }
        if(!incR.empty()){
            std::sort(incR.begin(), incR.end());
            double medR=incR[incR.size()/2];
            bgSize[vi]=std::min(bgSize[vi], medR);
        }
    }

    // Estimate triangle count using two components:
    //
    // 1) Background mesh: integrate 1/(alpha*H^2) over domain using CDT triangles
    //    where H is the background (area-constraint) mesh size.
    //
    // 2) Refinement seeds: vertices where lfs << H are local refinement points
    //    (intentionally close vertices to force mesh density at corners, etc.)
    //    Model each seed as a circle of influence with radius r = c*H where
    //    mesh grades from h=lfs at center to h=H at radius r.
    //    The grading follows h(d) = lfs + d*(H-lfs)/r for distance d from seed.
    //    Integrating 1/h^2 over this circle (in polar coords):
    //      contribution = integral_0^r (2*pi*d / (alpha*h(d)^2)) dd
    //    Approximation: the dense core (radius ~lfs) contributes pi*lfs^2/(alpha*lfs^2)
    //    = pi/alpha ≈ 7.3 triangles, plus the grading annulus contributes
    //    ~2*pi/alpha * ln(H/lfs) triangles (from the 1/h^2 integral with linear grading).
    //    Total per seed ≈ (pi/alpha) * (1 + 2*ln(H/lfs))
    //    Subtract the background contribution for the circle area (pi*c^2*H^2/(alpha*H^2)
    //    = pi*c^2/alpha) to avoid double-counting.

    double alpha=std::sqrt(3.0)/4.0;
    double piOverAlpha=M_PI/alpha; // ≈ 7.26
    double cosMinAngle=std::cos(minAng*M_PI/180.0);

    // Pass 1: background mesh estimate from CDT triangles using H
    double estBackground=0;
    int nCDT=0, nBadCDT=0;
    double totalArea=0;
    for(int i=0;i<nt;i++){
        auto& t=mesh.triangles[i];
        if(t.v[0]<0) continue;
        nCDT++;
        double area=triArea(mesh.vertices[t.v[0]],mesh.vertices[t.v[1]],mesh.vertices[t.v[2]]);
        totalArea+=area;
        double H=(bgSize[t.v[0]]+bgSize[t.v[1]]+bgSize[t.v[2]])/3.0;
        double bgEst=area/(alpha*H*H);
        // Check for bad angle
        auto& p0=mesh.vertices[t.v[0]];
        auto& p1=mesh.vertices[t.v[1]];
        auto& p2=mesh.vertices[t.v[2]];
        double e0=std::sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
        double e1=std::sqrt((p0.x-p2.x)*(p0.x-p2.x)+(p0.y-p2.y)*(p0.y-p2.y));
        double e2=std::sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y));
        double cosA0=(e1*e1+e2*e2-e0*e0)/(2*e1*e2+1e-30);
        double cosA1=(e0*e0+e2*e2-e1*e1)/(2*e0*e2+1e-30);
        double cosA2=(e0*e0+e1*e1-e2*e2)/(2*e0*e1+1e-30);
        double cosMax=std::max({cosA0,cosA1,cosA2});
        bool isBad=(cosMax>cosMinAngle);
        if(isBad) nBadCDT++;
        // When area constraints force finer mesh AND the triangle has bad angles,
        // the area-sized output triangles also need angle refinement, roughly
        // doubling the count vs area-only.  Apply 1.7x interaction factor.
        if(isBad && bgEst>6.0)
            estBackground+=bgEst*1.7;
        else
            estBackground+=std::max(bgEst, 6.0); // min 6 per bad CDT tri
    }

    // Pass 2: refinement seed contributions
    double estSeeds=0;
    int nSeeds=0;
    for(int vi=0;vi<nv;vi++){
        double h=lfs[vi], H=bgSize[vi];
        if(h>=1e29 || h>=H*0.5) continue; // not a refinement seed
        nSeeds++;
        // Circle of influence: grading from h to H
        double ratio=H/h;
        double logR=std::log(ratio);
        // Triangles from grading: pi/alpha * (1 + 2*ln(ratio))
        // Correction: Ruppert's geometric grading creates quadratically more
        // triangles than the linear grading model at high H/h ratios.
        double seedContrib=piOverAlpha*(1.0+2.0*logR)*std::max(1.0, logR);
        // Subtract background that would have been there anyway
        // Background in circle area ~pi*(c*H)^2 / (alpha*H^2) = pi*c^2/alpha
        // With c~1 (grading distance ≈ H), subtract pi/alpha ≈ 7.3
        seedContrib-=piOverAlpha;
        if(seedContrib>0) estSeeds+=seedContrib;
    }

    // Pass 3: segment LFS contributions
    // A segment with LFS constraint h generates a grading band on both sides.
    // Theoretical integral (linear grading) telescopes to 2*L/(alpha*h), but
    // Ruppert's geometric grading is ~3x denser, giving 6*L/(alpha*h).
    // NOTE: this assumes grading on both sides (interior segment).  Boundary
    // segments only grade on one side, so this overestimates by ~2x for them.
    // Conservative is fine for a safety cap; if a tighter estimate is ever
    // needed, halve the contribution for boundary segments (those with only
    // one incident triangle on the constrained edge).
    double estSegLfs=0;
    for(auto& s:mesh.segments){
        if(s.lfs<=0) continue;
        double Ldx=mesh.vertices[s.v1].x-mesh.vertices[s.v0].x, Ldy=mesh.vertices[s.v1].y-mesh.vertices[s.v0].y;
        double L=std::sqrt(Ldx*Ldx+Ldy*Ldy);
        estSegLfs+=6.0*L/(alpha*s.lfs);
    }

    double estTotalTris=estBackground+estSeeds+estSegLfs;
    // Steiner limit: 5x estimate gives generous headroom for legitimate cases
    // while stopping pathological non-convergence
    int maxSteiner=std::max(10000, (int)(estTotalTris*5.0));
    if(!opts.quiet){
        std::cerr<<"  CDT: "<<nCDT<<" tris ("<<nBadCDT<<" bad), area="
                 <<totalArea<<", seeds: "<<nSeeds<<"\n";
        std::cerr<<"  Est triangles: "<<(int)estTotalTris
                 <<", Steiner limit: "<<maxSteiner<<"\n";
    }
    int steinerCount=0;

    EdgeSet constrainedEdges;
    for(auto& s:mesh.segments) constrainedEdges.insert(edgeKey(s.v0,s.v1));
    double gxmin=1e30,gxmax=-1e30,gymin=1e30,gymax=-1e30;
    for(auto& v:mesh.vertices){
        gxmin=std::min(gxmin,v.x); gxmax=std::max(gxmax,v.x);
        gymin=std::min(gymin,v.y); gymax=std::max(gymax,v.y);
    }
    gxmin-=1; gymin-=1; gxmax+=1; gymax+=1;
    int gNX=50, gNY=50;
    double gcell=std::max((gxmax-gxmin)/gNX,(gymax-gymin)/gNY);
    gNX=(int)((gxmax-gxmin)/gcell)+1;
    gNY=(int)((gymax-gymin)/gcell)+1;

    std::unordered_map<int,std::vector<int>> segGrid;

    auto gridAddSeg=[&](int si){
        auto& s=mesh.segments[si];
        auto& sa=mesh.vertices[s.v0]; auto& sb=mesh.vertices[s.v1];
        double mx=(sa.x+sb.x)*0.5, my=(sa.y+sb.y)*0.5;
        double r2=(sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y);
        double r=std::sqrt(r2)*0.5;
        int cx0=std::max(0,(int)((mx-r-gxmin)/gcell));
        int cy0=std::max(0,(int)((my-r-gymin)/gcell));
        int cx1=std::min(gNX-1,(int)((mx+r-gxmin)/gcell));
        int cy1=std::min(gNY-1,(int)((my+r-gymin)/gcell));
        for(int cx=cx0;cx<=cx1;cx++)
            for(int cy=cy0;cy<=cy1;cy++)
                segGrid[cx*gNY+cy].push_back(si);
    };
    auto gridRemoveSeg=[&](int si){
        auto& s=mesh.segments[si];
        auto& sa=mesh.vertices[s.v0]; auto& sb=mesh.vertices[s.v1];
        double mx=(sa.x+sb.x)*0.5, my=(sa.y+sb.y)*0.5;
        double r2=(sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y);
        double r=std::sqrt(r2)*0.5;
        int cx0=std::max(0,(int)((mx-r-gxmin)/gcell));
        int cy0=std::max(0,(int)((my-r-gymin)/gcell));
        int cx1=std::min(gNX-1,(int)((mx+r-gxmin)/gcell));
        int cy1=std::min(gNY-1,(int)((my+r-gymin)/gcell));
        for(int cx=cx0;cx<=cx1;cx++)
            for(int cy=cy0;cy<=cy1;cy++){
                auto& v=segGrid[cx*gNY+cy];
                v.erase(std::remove(v.begin(),v.end(),si),v.end());
            }
    };

    for(int si=0;si<(int)mesh.segments.size();si++) gridAddSeg(si);

    // Map from edge key to segment index, for finding PBC partner segments
    std::map<std::pair<int,int>,int> edgeToSeg;
    for(int si=0;si<(int)mesh.segments.size();si++)
        edgeToSeg[edgeKey(mesh.segments[si].v0,mesh.segments[si].v1)]=si;

    // Core segment split helper (no PBC synchronization)
    auto splitSegmentCore=[&](int si, std::vector<int>& outAffected, int hintTri=-1) -> int {
        auto seg=mesh.segments[si]; // copy before modifying
        auto& sa=mesh.vertices[seg.v0]; auto& sb=mesh.vertices[seg.v1];
        double mx=(sa.x+sb.x)*0.5, my=(sa.y+sb.y)*0.5;
        double midLfs = seg.lfs;
        Point mid{mx,my,(int)mesh.vertices.size(),seg.marker,{},midLfs};

        auto ek=edgeKey(seg.v0,seg.v1);
        constrainedEdges.erase(ek);
        edgeToSeg.erase(ek);
        gridRemoveSeg(si);

        int midIdx=insertPointLawson(mesh,mid,hintTri,constrainedEdges,outAffected);
        if(midIdx<0){
            constrainedEdges.insert(ek);
            edgeToSeg[ek]=si;
            gridAddSeg(si);
            return -1;
        }
        steinerCount++;
        Segment s1{seg.v0,midIdx,seg.marker,seg.lfs,seg.pbc_partner,seg.pbc_type};
        Segment s2{midIdx,seg.v1,seg.marker,seg.lfs,seg.pbc_partner,seg.pbc_type};
        mesh.segments[si]=s1;
        int newSi=(int)mesh.segments.size();
        mesh.segments.push_back(s2);
        constrainedEdges.insert(edgeKey(s1.v0,s1.v1));
        constrainedEdges.insert(edgeKey(s2.v0,s2.v1));
        edgeToSeg[edgeKey(s1.v0,s1.v1)]=si;
        edgeToSeg[edgeKey(s2.v0,s2.v1)]=newSi;
        gridAddSeg(si);
        gridAddSeg(newSi);
        return midIdx;
    };

    // Lambda to split a segment at its midpoint using Lawson insertion.
    // For PBC segments, also split the partner segment to keep boundaries
    // synchronized ("tube" topology — both boundaries are conceptually the
    // same edge, so splitting one splits the other).
    auto splitSegment=[&](int si, std::vector<int>& outAffected, int hintTri=-1) -> bool {
        auto seg=mesh.segments[si]; // copy — seg may be modified by splitSegmentCore
        int midIdx=splitSegmentCore(si, outAffected, hintTri);
        if(midIdx<0) return false;

        // Synchronized PBC partner split: if this is a PBC segment, find and
        // split the partner. Suppress cascade by not adding partner's affected
        // triangles to the quality check.
        if(seg.pbc_type>=0){
            int tv0=mesh.pbc_twin.count(seg.v0)?mesh.pbc_twin[seg.v0]:-1;
            int tv1=mesh.pbc_twin.count(seg.v1)?mesh.pbc_twin[seg.v1]:-1;
            if(tv0>=0 && tv1>=0){
                auto pek=edgeKey(tv0,tv1);
                auto pit=edgeToSeg.find(pek);
                if(pit!=edgeToSeg.end()){
                    int partnerSi=pit->second;
                    std::vector<int> partnerAffected;
                    int partnerMid=splitSegmentCore(partnerSi, partnerAffected);
                    if(partnerMid>=0){
                        // Link the two midpoints as PBC twins
                        mesh.pbc_twin[midIdx]=partnerMid;
                        mesh.pbc_twin[partnerMid]=midIdx;
                        mesh.pbc_node_type[midIdx]=seg.pbc_type;
                        mesh.pbc_node_type[partnerMid]=seg.pbc_type;
                        // DO NOT add partnerAffected triangles to the quality PQ.
                        // This breaks the encroachment cascade. Any bad triangles
                        // near the partner will be caught in subsequent iterations.
                    }
                }
            }
        }
        return true;
    };

    // Check if point (px,py) encroaches any segment; return segment index or -1
    auto findEncroached=[&](double px, double py) -> int {
        int gcx=std::max(0,std::min(gNX-1,(int)((px-gxmin)/gcell)));
        int gcy=std::max(0,std::min(gNY-1,(int)((py-gymin)/gcell)));
        // Check nearby grid cells
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++){
            int cx=gcx+dx, cy=gcy+dy;
            if(cx<0||cx>=gNX||cy<0||cy>=gNY) continue;
            auto it=segGrid.find(cx*gNY+cy);
            if(it==segGrid.end()) continue;
            for(int si:it->second){
                auto& s=mesh.segments[si];
                auto& sa=mesh.vertices[s.v0]; auto& sb=mesh.vertices[s.v1];
                double mx2=(sa.x+sb.x)*0.5, my2=(sa.y+sb.y)*0.5;
                double sr2=((sb.x-sa.x)*(sb.x-sa.x)+(sb.y-sa.y)*(sb.y-sa.y))*0.25;
                double dd2=(px-mx2)*(px-mx2)+(py-my2)*(py-my2);
                if(dd2<sr2) return si;
            }
        }
        return -1;
    };

    using PQItem = std::tuple<double, int, int>; // metric, triIdx, generation
    std::priority_queue<PQItem> pq;

    for(int i=0;i<(int)mesh.triangles.size();i++){
        auto mr = triMetric(mesh, i, minAng, globalMaxA, useRegionArea, cosMinAng);
        if(mr.metric > 0) pq.push({mr.metric, i, mesh.triangles[i].generation});
    }

    // Precompute off-center constant (loop-invariant)
    const double offconst=minAng>0 ? 0.475*std::sqrt((1.0+std::cos(minAng*M_PI/180.0))/(1.0-std::cos(minAng*M_PI/180.0))) : 0;

    // Pre-allocated buffers reused across iterations
    std::vector<int> affected, toCheck, segAff;
    std::vector<std::pair<int,int>> flipBuf;

    int dbgEncroach=0, dbgInsert=0, dbgSkip=0, dbgFail=0;
    while(!pq.empty() && steinerCount<maxSteiner){
        auto [metric, badIdx, gen] = pq.top();
        pq.pop();

        if(badIdx >= (int)mesh.triangles.size()) continue;
        auto& bt = mesh.triangles[badIdx];
        if(bt.v[0]<0) continue;
        // Fast staleness check: skip if triangle was modified since this PQ entry
        if(bt.generation != gen) continue;
        auto tmr = triMetric(mesh, badIdx, minAng, globalMaxA, useRegionArea, cosMinAng);
        if(tmr.metric <= 0) continue;
        bool angViol=tmr.angViol, areaViol=tmr.areaViol;
        double area=tmr.area, effectiveMaxArea=tmr.effectiveMaxArea;

        auto& a=mesh.vertices[bt.v[0]];
        auto& b=mesh.vertices[bt.v[1]];
        auto& c=mesh.vertices[bt.v[2]];

        // Skip triangles where both edges at the smallest-angle vertex are segments
        // (the small angle is a segment intersection that can't be fixed).
        // This matches Triangle 1.3's testtriangle() logic.
        if(angViol && !areaViol){
            // Find vertex with smallest angle = vertex opposite the shortest edge
            // Use cached squared edge lengths from triMetric
            double la2=tmr.la2, lb2=tmr.lb2, lc2=tmr.lc2;
            int minV=0; // opposite la (edge b-c)
            if(lb2<la2 && lb2<lc2) minV=1; // opposite lb (edge a-c)
            else if(lc2<la2) minV=2;        // opposite lc (edge a-b)
            // Check if both edges at the smallest-angle vertex are segments
            int va=bt.v[minV], vb=bt.v[(minV+1)%3], vc=bt.v[(minV+2)%3];
            if(constrainedEdges.count(edgeKey(va,vb)) && constrainedEdges.count(edgeKey(va,vc))){
                dbgSkip++; continue;
            }
        }

        // Compute circumcenter or off-center
        double tcx,tcy;
        if(angViol){
            // Inline circumcenter using cached values from triMetric
            double xdo=b.x-a.x, ydo=b.y-a.y; // a->b
            double xao=c.x-a.x, yao=c.y-a.y; // a->c
            double D=2.0*(xdo*yao-ydo*xao);   // = 2*cross from triMetric
            if(std::abs(D)<EPS){ tcx=(a.x+b.x+c.x)/3.0; tcy=(a.y+b.y+c.y)/3.0; }
            else {
                // dodist=xdo²+ydo²=lc2, aodist=xao²+yao²=lb2 (from triMetric)
                double dodist=tmr.lc2, aodist=tmr.lb2;
                tcx=a.x+(yao*dodist-ydo*aodist)/D;
                tcy=a.y+(xdo*aodist-xao*dodist)/D;
            }
            // Off-center: limit distance from shortest edge (Üngör)
            double dodist=tmr.lc2, aodist=tmr.lb2, dadist=tmr.la2;
            double dx=tcx-a.x, dy=tcy-a.y; // circumcenter relative to a
            if(dodist<aodist && dodist<dadist){
                // Shortest edge is a-b (origin-dest)
                double dxoff=0.5*xdo - offconst*ydo;
                double dyoff=0.5*ydo + offconst*xdo;
                if(dxoff*dxoff+dyoff*dyoff < dx*dx+dy*dy){
                    tcx=a.x+dxoff; tcy=a.y+dyoff;
                }
            } else if(aodist<dadist){
                // Shortest edge is a-c (origin-apex)
                double dxoff=0.5*xao + offconst*yao;
                double dyoff=0.5*yao - offconst*xao;
                if(dxoff*dxoff+dyoff*dyoff < dx*dx+dy*dy){
                    tcx=a.x+dxoff; tcy=a.y+dyoff;
                }
            } else {
                // Shortest edge is b-c (dest-apex)
                double dxoff=0.5*(c.x-b.x) - offconst*(c.y-b.y);
                double dyoff=0.5*(c.y-b.y) + offconst*(c.x-b.x);
                if(dxoff*dxoff+dyoff*dyoff < (dx-xdo)*(dx-xdo)+(dy-ydo)*(dy-ydo)){
                    tcx=b.x+dxoff; tcy=b.y+dyoff;
                }
            }
        } else {
            tcx=(a.x+b.x+c.x)/3.0;
            tcy=(a.y+b.y+c.y)/3.0;
        }

        // Pre-insertion encroachment check: if circumcenter/off-center encroaches
        // any segment, split that segment instead of inserting the point.
        // With -Y, skip if the encroached segment is a boundary segment (marker!=0).
        int encSeg=findEncroached(tcx,tcy);
        if(encSeg>=0 && opts.no_steiner && mesh.segments[encSeg].marker!=0){
            continue; // -Y: don't add Steiner points on boundary segments
        }
        if(encSeg>=0){
            segAff.clear();
            if(splitSegment(encSeg, segAff, badIdx)){
                dbgEncroach++;
                {auto mr=triMetric(mesh,badIdx,minAng,globalMaxA,useRegionArea,cosMinAng);
                if(mr.metric>0) pq.push({mr.metric,badIdx,mesh.triangles[badIdx].generation});}
                for(int ti2:segAff){
                    if(ti2>=0&&ti2<(int)mesh.triangles.size()&&mesh.triangles[ti2].v[0]>=0)
                        for(int k=0;k<3;k++){
                            int nb=mesh.triangles[ti2].neighbors[k];
                            if(nb>=0){ auto mr2=triMetric(mesh,nb,minAng,globalMaxA,useRegionArea,cosMinAng); if(mr2.metric>0) pq.push({mr2.metric,nb,mesh.triangles[nb].generation}); }
                        }
                }
            }
            continue;
        }

        // Insert off-center/circumcenter using Lawson flips
        Point np{tcx,tcy,(int)mesh.vertices.size(),0,{}};
        affected.clear();
        int newPt=insertPointLawson(mesh,np,badIdx,constrainedEdges,affected,&flipBuf);
        if(newPt<0){ dbgFail++; continue; }

        steinerCount++;
        dbgInsert++;

        // Re-check affected triangles and their neighbors
        // Use a flat sorted vector instead of unordered_set to avoid heap overhead
        toCheck.assign(affected.begin(), affected.end());
        for(int ti : affected){
            if(ti>=0 && ti<(int)mesh.triangles.size() && mesh.triangles[ti].v[0]>=0){
                for(int j=0;j<3;j++){
                    int nb=mesh.triangles[ti].neighbors[j];
                    if(nb>=0) toCheck.push_back(nb);
                }
            }
        }
        std::sort(toCheck.begin(), toCheck.end());
        toCheck.erase(std::unique(toCheck.begin(), toCheck.end()), toCheck.end());
        for(int ti : toCheck){
            if(ti<0 || ti>=(int)mesh.triangles.size()) continue;
            auto mr=triMetric(mesh,ti,minAng,globalMaxA,useRegionArea,cosMinAng);
            if(mr.metric>0) pq.push({mr.metric,ti,mesh.triangles[ti].generation});
        }
    }
    if(!opts.quiet)
        std::cerr<<"Quality refinement: "<<dbgInsert<<" Steiner points, "
                 <<dbgEncroach<<" segment splits, "<<dbgSkip<<" skipped, "
                 <<dbgFail<<" failed insertions\n";

}

// ============================================================
// Edge extraction
// ============================================================

void extractEdges(Mesh& mesh){
    mesh.edges.clear();
    EdgeSet seen;
    for(auto& t:mesh.triangles){
        if(t.v[0]<0) continue;
        for(int j=0;j<3;j++){
            auto key=edgeKey(t.v[j],t.v[(j+1)%3]);
            if(!seen.count(key)){seen.insert(key); mesh.edges.push_back({t.v[j],t.v[(j+1)%3]});}
        }
    }
}

// ============================================================
// File I/O
// ============================================================

static void skipLine(std::istream& in){
    in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}
static void skipComments(std::istream& in){
    // Skip whitespace, then skip lines starting with '#'
    in>>std::ws;
    while(in.peek()=='#'){ skipLine(in); in>>std::ws; }
}

bool readNodeFile(const std::string& filename, std::vector<Point>& pts, int& nAttribs, int& nMarkers){
    std::ifstream f(filename);
    if(!f){std::cerr<<"Cannot open "<<filename<<"\n";return false;}
    skipComments(f);
    int n,dim; f>>n>>dim>>nAttribs>>nMarkers;
    pts.resize(n);
    for(int i=0;i<n;i++){
        skipComments(f);
        int idx; f>>idx>>pts[i].x>>pts[i].y;
        pts[i].id=idx;
        pts[i].attribs.resize(nAttribs);
        for(int a=0;a<nAttribs;a++) f>>pts[i].attribs[a];
        if(nMarkers>0) f>>pts[i].marker; else pts[i].marker=0;
        // Optional LFS column after marker
        std::string rest; std::getline(f,rest);
        std::istringstream rss(rest);
        double lval;
        if(rss>>lval) pts[i].lfs=lval;
    }
    return true;
}

bool readPolyFile(const std::string& filename, Mesh& mesh, int& nAttribs){
    auto& pts=mesh.vertices;
    auto& segs=mesh.segments;
    auto& holes=mesh.holes;
    auto& regions=mesh.regions;

    std::ifstream f(filename);
    if(!f){std::cerr<<"Cannot open "<<filename<<"\n";return false;}
    skipComments(f);
    int nv,dim,nmark; f>>nv>>dim>>nAttribs>>nmark;
    pts.resize(nv);
    int firstVertIdx=std::numeric_limits<int>::max();
    for(int i=0;i<nv;i++){
        skipComments(f);
        int idx; f>>idx>>pts[i].x>>pts[i].y;
        pts[i].id=idx;
        firstVertIdx=std::min(firstVertIdx,idx);
        pts[i].attribs.resize(nAttribs);
        for(int a=0;a<nAttribs;a++) f>>pts[i].attribs[a];
        if(nmark>0) f>>pts[i].marker; else pts[i].marker=0;
        // Optional LFS column after marker
        std::string rest; std::getline(f,rest);
        std::istringstream rss(rest);
        double lval;
        if(rss>>lval) pts[i].lfs=lval;
    }
    int base=(firstVertIdx==0)?0:1;

    skipComments(f);
    int ns,smark; f>>ns>>smark;
    segs.resize(ns);
    for(int i=0;i<ns;i++){
        skipComments(f);
        int idx; f>>idx>>segs[i].v0>>segs[i].v1;
        segs[i].v0-=base; segs[i].v1-=base;
        if(smark>0) f>>segs[i].marker; else segs[i].marker=0;
        // Optional LFS column after marker
        std::string rest; std::getline(f,rest);
        std::istringstream rss(rest);
        double lval;
        if(rss>>lval) segs[i].lfs=lval;
    }
    skipComments(f);
    int nh; f>>nh;
    holes.resize(nh);
    for(int i=0;i<nh;i++){
        skipComments(f);
        int idx; f>>idx>>holes[i].x>>holes[i].y;
    }
    skipComments(f);
    {
        int nr=0;
        if(f>>nr){
            regions.resize(nr);
            for(int i=0;i<nr;i++){
                skipComments(f);
                int idx; f>>idx>>regions[i].x>>regions[i].y>>regions[i].attrib;
                regions[i].max_area=-1.0;
                std::string rest; std::getline(f,rest);
                std::istringstream rss(rest);
                double aval;
                if(rss>>aval) regions[i].max_area=aval;
            }
        }
    }
    // Clear stream state so optional extension sections can be attempted
    f.clear();

    // ---- tangle extensions: arcs and PBCs (optional, after regions) ----

    // Arc segments
    {
        skipComments(f);
        int nArcs=0, arcHasMarkers=0;
        if(f>>nArcs>>arcHasMarkers){
            for(int i=0;i<nArcs;i++){
                skipComments(f);
                int idx, v0, v1;
                double angle, maxSegAngle;
                f>>idx>>v0>>v1>>angle>>maxSegAngle;
                v0-=base; v1-=base;
                int marker=0;
                if(arcHasMarkers) f>>marker;
                // Optional LFS after marker
                double arcLfs=-1.0;
                std::string rest; std::getline(f,rest);
                std::istringstream rss(rest);
                double lval;
                if(rss>>lval) arcLfs=lval;

                // Sanity checks
                angle=std::abs(angle);
                if(angle<1.0||angle>180.0){
                    std::cerr<<"Warning: arc "<<idx<<" angle "<<angle
                             <<"° clamped to [1,180]\n";
                    angle=std::max(1.0, std::min(180.0, angle));
                }
                if(maxSegAngle>10.0||maxSegAngle<=0.01) maxSegAngle=10.0;

                // Discretize arc into chord segments
                double ax=pts[v0].x, ay=pts[v0].y;
                double bx=pts[v1].x, by=pts[v1].y;
                double d=std::sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
                if(d<1e-30){ segs.push_back({v0, v1, marker, arcLfs}); continue; }

                double tta=angle*M_PI/180.0;
                double R=d/(2.0*std::sin(tta/2.0));
                double tx=(bx-ax)/d, ty=(by-ay)/d;
                double h=std::sqrt(std::max(0.0, R*R - d*d/4.0));
                double cx=ax + (d/2.0)*tx + h*(-ty);
                double cy=ay + (d/2.0)*ty + h*(tx);

                int k=(int)std::ceil(angle/maxSegAngle);
                if(k<1) k=1;

                if(k==1){
                    segs.push_back({v0, v1, marker, arcLfs});
                } else {
                    double dtheta=angle*M_PI/(180.0*(double)k);
                    double cosD=std::cos(dtheta), sinD=std::sin(dtheta);
                    double px=ax, py=ay;
                    int prevIdx=v0;
                    for(int j=0;j<k;j++){
                        double dx2=px-cx, dy2=py-cy;
                        px=cx + dx2*cosD - dy2*sinD;
                        py=cy + dx2*sinD + dy2*cosD;
                        if(j<k-1){
                            int newIdx=(int)pts.size();
                            pts.push_back({px, py, newIdx, marker, {}, arcLfs});
                            segs.push_back({prevIdx, newIdx, marker, arcLfs});
                            prevIdx=newIdx;
                        } else {
                            segs.push_back({prevIdx, v1, marker, arcLfs});
                        }
                    }
                }
            }
        }
    }
    f.clear();

    // PBC definitions: pair two boundary markers
    {
        skipComments(f);
        int nPBCs=0;
        if(f>>nPBCs){
            for(int i=0;i<nPBCs;i++){
                skipComments(f);
                int idx, markerA, markerB, type;
                f>>idx>>markerA>>markerB>>type;
                // Mark all segments with markerA or markerB as PBC
                for(auto& s:segs){
                    if(s.marker==markerA || s.marker==markerB){
                        s.pbc_type=type;
                    }
                }
            }
        }
    }

    return true;
}



void writeNodeFile(const std::string& fn, const std::vector<Point>& pts,
                   int nAttribs, bool hasMarkers, const Options& opts){
    std::ofstream f(fn);
    f<<std::setprecision(17);
    bool hasLfs=false;
    for(auto& p:pts) if(p.lfs>0){hasLfs=true;break;}
    f<<pts.size()<<" 2 "<<nAttribs<<" "<<(hasMarkers?1:0)<<"\n";
    for(int i=0;i<(int)pts.size();i++){
        int outIdx=opts.zero_indexed?i:(i+1);
        f<<outIdx<<" "<<pts[i].x<<" "<<pts[i].y;
        for(double a:pts[i].attribs) f<<" "<<a;
        if(hasMarkers) f<<" "<<pts[i].marker;
        if(hasLfs) f<<" "<<pts[i].lfs;
        f<<"\n";
    }
}

void writeEleFile(const std::string& fn, const Mesh& mesh, bool hasAttrib, const Options& opts){
    std::ofstream f(fn);
    f<<mesh.triangles.size()<<" 3 "<<(hasAttrib?1:0)<<"\n";
    for(int i=0;i<(int)mesh.triangles.size();i++){
        auto& t=mesh.triangles[i];
        int outIdx=opts.zero_indexed?i:(i+1);
        int v0=opts.zero_indexed?t.v[0]:(t.v[0]+1);
        int v1=opts.zero_indexed?t.v[1]:(t.v[1]+1);
        int v2=opts.zero_indexed?t.v[2]:(t.v[2]+1);
        f<<outIdx<<" "<<v0<<" "<<v1<<" "<<v2;
        if(hasAttrib) f<<" "<<t.region_attrib;
        f<<"\n";
    }
}

void writeEdgeFile(const std::string& fn, const Mesh& mesh, const Options& opts){
    std::ofstream f(fn);
    std::map<std::pair<int,int>,int> edgeCount;
    for(auto& t:mesh.triangles){
        if(t.v[0]<0) continue;
        for(int j=0;j<3;j++) edgeCount[edgeKey(t.v[j],t.v[(j+1)%3])]++;
    }
    // Build map from edge to segment marker
    std::map<std::pair<int,int>,int> segMarker;
    for(auto& s:mesh.segments) segMarker[edgeKey(s.v0,s.v1)]=s.marker;
    f<<mesh.edges.size()<<" 1\n";
    for(int i=0;i<(int)mesh.edges.size();i++){
        int a=mesh.edges[i].first, b=mesh.edges[i].second;
        int outIdx=opts.zero_indexed?i:(i+1);
        int oa=opts.zero_indexed?a:(a+1);
        int ob=opts.zero_indexed?b:(b+1);
        auto key=edgeKey(a,b);
        int marker=0;
        bool isBoundary = (edgeCount[key]==1);
        auto it=segMarker.find(key);
        if(it!=segMarker.end()){
            marker=it->second;
            if(marker==0 && isBoundary) marker=1; // boundary segments get at least 1
        } else if(isBoundary){
            marker=1;
        }
        f<<outIdx<<" "<<oa<<" "<<ob<<" "<<marker<<"\n";
    }
}

void writeNeighFile(const std::string& fn, const Mesh& mesh, const Options& opts){
    std::ofstream f(fn);
    f<<mesh.triangles.size()<<" 3\n";
    for(int i=0;i<(int)mesh.triangles.size();i++){
        auto& t=mesh.triangles[i];
        int outIdx=opts.zero_indexed?i:(i+1);
        auto nb=[&](int n)->int{return(n<0)?-1:(opts.zero_indexed?n:(n+1));};
        f<<outIdx<<" "<<nb(t.neighbors[0])<<" "<<nb(t.neighbors[1])<<" "<<nb(t.neighbors[2])<<"\n";
    }
}

void writePolyFile(const std::string& fn, const Mesh& mesh, const Options& opts){
    std::ofstream f(fn);
    f<<std::setprecision(17);
    bool hasLfs=false;
    for(auto& p:mesh.vertices) if(p.lfs>0){hasLfs=true;break;}
    bool hasSegLfs=false;
    for(auto& s:mesh.segments) if(s.lfs>0){hasSegLfs=true;break;}
    int nv=(int)mesh.vertices.size();
    f<<nv<<" 2 0 1\n";
    for(int i=0;i<nv;i++){
        int outIdx=opts.zero_indexed?i:(i+1);
        f<<outIdx<<" "<<mesh.vertices[i].x<<" "<<mesh.vertices[i].y<<" "<<mesh.vertices[i].marker;
        if(hasLfs) f<<" "<<mesh.vertices[i].lfs;
        f<<"\n";
    }
    f<<mesh.segments.size()<<" 1\n";
    for(int i=0;i<(int)mesh.segments.size();i++){
        int outIdx=opts.zero_indexed?i:(i+1);
        int oa=opts.zero_indexed?mesh.segments[i].v0:(mesh.segments[i].v0+1);
        int ob=opts.zero_indexed?mesh.segments[i].v1:(mesh.segments[i].v1+1);
        f<<outIdx<<" "<<oa<<" "<<ob<<" "<<mesh.segments[i].marker;
        if(hasSegLfs) f<<" "<<mesh.segments[i].lfs;
        f<<"\n";
    }
    f<<"0\n";
}

void writePbcFile(const std::string& fn, const Mesh& mesh, const Options& opts){
    std::ofstream f(fn);
    f<<std::setprecision(17);
    f<<mesh.pbc_pairs.size()<<"\n";
    for(int i=0;i<(int)mesh.pbc_pairs.size();i++){
        auto& p=mesh.pbc_pairs[i];
        f<<i<<" "<<p.node_a<<" "<<p.node_b<<" "<<p.type<<"\n";
    }
}

// ============================================================
// Argument parsing
// ============================================================

Options parseOptions(const std::string& switchStr){
    Options opts;
    int i=0;
    while(i<(int)switchStr.size()){
        char c=switchStr[i++];
        switch(c){
            case 'p': opts.pslg=true; break;
            case 'P': opts.no_poly_out=true; break;
            case 'j': opts.jettison=true; break;
            case 'q':
                opts.quality=true;
                if(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.')){
                    std::string num;
                    while(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.'))
                        num+=switchStr[i++];
                    opts.min_angle=std::stod(num);
                } break;
            case 'e': opts.edges=true; break;
            case 'A': opts.regions=true; break;
            case 'a':
                opts.area_limit=true;
                if(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.')){
                    std::string num;
                    while(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.'))
                        num+=switchStr[i++];
                    opts.max_area=std::stod(num);
                } break;
            case 'z': opts.zero_indexed=true; opts.first_index=0; break;
            case 'Q': opts.quiet=true; break;
            case 'I': opts.suppress_iter=true; break;
            case 'O': opts.no_holes=true; break;
            case 'Y': opts.no_steiner=true; break;
            case 'n': opts.neighbors=true; break;
            case 'C':
                opts.clean_pslg=true;
                if(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.')){
                    std::string num;
                    while(i<(int)switchStr.size()&&(std::isdigit(switchStr[i])||switchStr[i]=='.'))
                        num+=switchStr[i++];
                    opts.clean_tol=std::stod(num);
                } break;
        }
    }
    return opts;
}

// ============================================================
// PSLG cleanup (-C flag): merge near-duplicate nodes, split
// segments at intersections, remove duplicate segments.
// Adapted from FEMM's EnforcePSLG().
// ============================================================

void cleanPSLG(Mesh& mesh, double tol, bool quiet){
    auto& pts = mesh.vertices;
    auto& segs = mesh.segments;
    int nv = (int)pts.size();

    // Compute auto tolerance if not specified
    if(tol < 0){
        double minX=pts[0].x, maxX=pts[0].x, minY=pts[0].y, maxY=pts[0].y;
        for(auto& p : pts){
            minX=std::min(minX,p.x); maxX=std::max(maxX,p.x);
            minY=std::min(minY,p.y); maxY=std::max(maxY,p.y);
        }
        tol = std::sqrt((maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY)) * 1e-6;
    }
    double tol2 = tol * tol;

    int mergedNodes=0, splitSegs=0, removedSegs=0;

    // 1. Merge near-duplicate nodes
    std::vector<int> remap(nv);
    std::iota(remap.begin(), remap.end(), 0);
    for(int i=0; i<nv; i++){
        if(remap[i] != i) continue; // already merged
        for(int j=i+1; j<nv; j++){
            if(remap[j] != j) continue;
            double dx=pts[j].x-pts[i].x, dy=pts[j].y-pts[i].y;
            if(dx*dx+dy*dy < tol2){
                remap[j] = i;
                mergedNodes++;
            }
        }
    }
    // Apply remap to segments
    for(auto& s : segs){
        s.v0 = remap[s.v0];
        s.v1 = remap[s.v1];
    }
    // Apply remap to holes and regions
    // (holes/regions use coordinates, not indices, so no remap needed)

    // 2. Remove degenerate segments (both endpoints same after merge)
    {
        std::vector<Segment> good;
        for(auto& s : segs)
            if(s.v0 != s.v1) good.push_back(s);
        removedSegs += (int)segs.size() - (int)good.size();
        segs = std::move(good);
    }

    // 3. Remove duplicate segments
    {
        std::set<std::pair<int,int>> seen;
        std::vector<Segment> good;
        for(auto& s : segs){
            auto k = edgeKey(s.v0, s.v1);
            if(seen.insert(k).second) good.push_back(s);
            else removedSegs++;
        }
        segs = std::move(good);
    }

    // 4. Split segments at near-coincident nodes (node lies on segment)
    {
        bool changed = true;
        while(changed){
            changed = false;
            for(int si=0; si<(int)segs.size(); si++){
                int v0=segs[si].v0, v1=segs[si].v1;
                double ax=pts[v0].x, ay=pts[v0].y;
                double bx=pts[v1].x, by=pts[v1].y;
                double dx=bx-ax, dy=by-ay;
                double len2=dx*dx+dy*dy;
                if(len2 < tol2) continue;
                for(int i=0; i<(int)pts.size(); i++){
                    if(i==v0 || i==v1) continue;
                    if(remap[i] != i && i < nv) continue; // merged away
                    double px=pts[i].x-ax, py=pts[i].y-ay;
                    // Distance from point to line segment
                    double t = (dx*px+dy*py)/len2;
                    if(t <= tol/std::sqrt(len2) || t >= 1.0-tol/std::sqrt(len2)) continue;
                    double cross = dx*py - dy*px;
                    if(cross*cross > tol2*len2) continue;
                    // Point i is on segment si — split it
                    Segment s2 = segs[si];
                    segs[si].v1 = i;
                    s2.v0 = i;
                    segs.push_back(s2);
                    splitSegs++;
                    changed = true;
                    break; // restart this segment (now shorter)
                }
            }
        }
    }

    // 5. Split segments at segment-segment intersections
    {
        bool changed = true;
        while(changed){
            changed = false;
            for(int i=0; i<(int)segs.size() && !changed; i++){
                double ax=pts[segs[i].v0].x, ay=pts[segs[i].v0].y;
                double bx=pts[segs[i].v1].x, by=pts[segs[i].v1].y;
                for(int j=i+1; j<(int)segs.size() && !changed; j++){
                    // Skip if they share an endpoint
                    if(segs[i].v0==segs[j].v0 || segs[i].v0==segs[j].v1 ||
                       segs[i].v1==segs[j].v0 || segs[i].v1==segs[j].v1) continue;
                    double cx=pts[segs[j].v0].x, cy=pts[segs[j].v0].y;
                    double dx2=pts[segs[j].v1].x, dy2=pts[segs[j].v1].y;
                    // Compute intersection
                    double d1x=bx-ax, d1y=by-ay;
                    double d2x=dx2-cx, d2y=dy2-cy;
                    double denom=d1x*d2y-d1y*d2x;
                    if(std::abs(denom) < 1e-15) continue;
                    double t=((cx-ax)*d2y-(cy-ay)*d2x)/denom;
                    double u=((cx-ax)*d1y-(cy-ay)*d1x)/denom;
                    if(t<=0 || t>=1 || u<=0 || u>=1) continue;
                    // Intersection at (ax+t*d1x, ay+t*d1y)
                    double ix=ax+t*d1x, iy=ay+t*d1y;
                    // Check if intersection is near an existing node
                    int nearNode=-1;
                    for(int k=0;k<(int)pts.size();k++){
                        double ex=pts[k].x-ix, ey=pts[k].y-iy;
                        if(ex*ex+ey*ey < tol2){ nearNode=k; break; }
                    }
                    if(nearNode<0){
                        nearNode=(int)pts.size();
                        pts.push_back({ix,iy,nearNode,0,{}});
                    }
                    // Split both segments at the intersection
                    Segment si2=segs[i]; si2.v0=nearNode;
                    segs[i].v1=nearNode;
                    segs.push_back(si2);
                    Segment sj2=segs[j]; sj2.v0=nearNode;
                    segs[j].v1=nearNode;
                    segs.push_back(sj2);
                    splitSegs+=2;
                    changed=true;
                }
            }
        }
    }

    // Compact: remove merged-away vertices and re-index
    if(mergedNodes > 0){
        std::vector<int> newIdx(pts.size(), -1);
        std::vector<Point> newPts;
        for(int i=0; i<(int)pts.size(); i++){
            int target = (i < nv) ? remap[i] : i;
            if(target != i) continue; // skip merged
            newIdx[i] = (int)newPts.size();
            pts[i].id = (int)newPts.size();
            newPts.push_back(pts[i]);
        }
        // Map merged indices through to their target's new index
        for(int i=0; i<nv; i++)
            if(remap[i] != i) newIdx[i] = newIdx[remap[i]];
        // Also map new intersection points
        for(int i=nv; i<(int)pts.size(); i++)
            if(newIdx[i]<0){ newIdx[i]=(int)newPts.size(); pts[i].id=(int)newPts.size(); newPts.push_back(pts[i]); }
        for(auto& s : segs){ s.v0=newIdx[s.v0]; s.v1=newIdx[s.v1]; }
        pts = std::move(newPts);
    }

    if(!quiet && (mergedNodes>0 || splitSegs>0 || removedSegs>0))
        std::cerr<<"  PSLG cleanup: merged "<<mergedNodes<<" nodes, split "
                 <<splitSegs<<" segments, removed "<<removedSegs<<" duplicates\n";
}

// ============================================================
// Main
// ============================================================

void printUsage(){
    std::cerr<<
"Usage: tangle [switches] inputfile\n"
"\n"
"Switches:\n"
"  -p    Read a Planar Straight-Line Graph (.poly file) and triangulate it.\n"
"  -P    Suppress output of .poly file.\n"
"  -j    Jettison vertices not present in any triangle from .node output.\n"
"  -q    Quality mesh generation (min angle, default 20 deg, e.g. -q28.5).\n"
"  -e    Output .edge file.\n"
"  -A    Assign regional attributes from .poly regions.\n"
"  -a    Area constraints (optional global max, e.g. -a0.5).\n"
"  -z    Number items starting from zero.\n"
"  -Q    Quiet mode.\n"
"  -I    Suppress iteration numbers on output file names.\n"
"  -Y    No Steiner points on boundary segments.\n"
"  -n    Output .neigh file.\n"
"  -C    Clean PSLG: merge near-duplicate nodes, split intersecting\n"
"        segments, remove duplicates. Optional tolerance (e.g. -C0.001);\n"
"        default is bounding-box diagonal * 1e-6.\n"
"\n"
"Input: .node or .poly files. Output: .1.node .1.ele [.1.edge] [.1.neigh] [.1.poly]\n";
}

int main(int argc, char* argv[]){
    if(argc<2){printUsage();return 1;}

    std::string switchStr, inputFile;
    for(int i=1;i<argc;i++){
        std::string arg=argv[i];
        if(arg[0]=='-') switchStr+=arg.substr(1);
        else inputFile=arg;
    }
    if(inputFile.empty()){std::cerr<<"No input file specified.\n";return 1;}

    Options opts=parseOptions(switchStr);

    std::string base=inputFile, ext;
    auto dot=base.rfind('.');
    if(dot!=std::string::npos){ext=base.substr(dot);base=base.substr(0,dot);}
    if(ext.empty()){ext=opts.pslg?".poly":".node";inputFile+=ext;}
    Mesh mesh;
    int nAttribs=0, nMarkers=0;

    if(opts.pslg||ext==".poly"){
        opts.pslg=true;
        if(!readPolyFile(inputFile,mesh,nAttribs))
            return 1;
        // Build PBC twin map for .poly PBC definitions.
        // PBC in .poly uses marker pairs: segments with marker_a pair with marker_b.
        // We find the two chains by marker and pair nodes by arc-length parameter.
        {
            // Collect distinct PBC marker groups
            std::set<int> pbcMarkers;
            for(auto& s:mesh.segments)
                if(s.pbc_type>=0) pbcMarkers.insert(s.marker);

            if(!pbcMarkers.empty()){
                auto buildChain=[&](const std::vector<int>& segIndices) -> std::vector<int> {
                    if(segIndices.empty()) return {};
                    std::map<int,std::vector<int>> nodeSegs;
                    for(int si:segIndices){
                        nodeSegs[mesh.segments[si].v0].push_back(si);
                        nodeSegs[mesh.segments[si].v1].push_back(si);
                    }
                    int startNode=-1;
                    for(auto& [node, slist]:nodeSegs)
                        if(slist.size()==1){ startNode=node; break; }
                    if(startNode<0) return {};
                    std::vector<int> chain; chain.push_back(startNode);
                    std::set<int> visited; int curNode=startNode;
                    while(true){
                        bool found=false;
                        for(int si:nodeSegs[curNode]){
                            if(visited.count(si)) continue;
                            visited.insert(si);
                            auto& s=mesh.segments[si];
                            int next=(s.v0==curNode)?s.v1:s.v0;
                            chain.push_back(next); curNode=next; found=true; break;
                        }
                        if(!found) break;
                    }
                    return chain;
                };

                // Group PBC markers into pairs by pbc_type
                // Each pair of markers with the same pbc_type forms a PBC boundary
                std::map<int,std::vector<int>> typeToMarkers; // pbc_type -> list of markers
                for(int m:pbcMarkers){
                    int pt=-1;
                    for(auto& s:mesh.segments)
                        if(s.marker==m && s.pbc_type>=0){ pt=s.pbc_type; break; }
                    if(pt>=0) typeToMarkers[pt].push_back(m);
                }

                for(auto& [pbcType, markers]:typeToMarkers){
                    // Process pairs of markers
                    for(int mi=0;mi+1<(int)markers.size();mi+=2){
                        int mA=markers[mi], mB=markers[mi+1];
                        std::vector<int> segsA, segsB;
                        for(int si=0;si<(int)mesh.segments.size();si++){
                            if(mesh.segments[si].marker==mA && mesh.segments[si].pbc_type>=0)
                                segsA.push_back(si);
                            if(mesh.segments[si].marker==mB && mesh.segments[si].pbc_type>=0)
                                segsB.push_back(si);
                        }
                        auto chainA=buildChain(segsA), chainB=buildChain(segsB);
                        if(chainA.empty()||chainB.empty()) continue;
                        if(chainA.size()!=chainB.size()){
                            std::cerr<<"Warning: PBC markers "<<mA<<" and "<<mB
                                     <<" have different chain lengths ("<<chainA.size()
                                     <<" vs "<<chainB.size()<<")\n";
                            continue;
                        }
                        // Orient chains consistently
                        auto dist2=[&](int a, int b){
                            double dx=mesh.vertices[a].x-mesh.vertices[b].x;
                            double dy=mesh.vertices[a].y-mesh.vertices[b].y;
                            return dx*dx+dy*dy;
                        };
                        double d_same=dist2(chainA[0],chainB[0])+dist2(chainA.back(),chainB.back());
                        double d_rev=dist2(chainA[0],chainB.back())+dist2(chainA.back(),chainB[0]);
                        if(d_rev<d_same) std::reverse(chainB.begin(), chainB.end());

                        for(int i=0;i<(int)chainA.size();i++){
                            int a=chainA[i], b=chainB[i];
                            if(a==b) continue;
                            mesh.pbc_twin[a]=b; mesh.pbc_twin[b]=a;
                            mesh.pbc_node_type[a]=pbcType;
                            mesh.pbc_node_type[b]=pbcType;
                        }
                    }
                }
            }
        }
    } else {
        if(!readNodeFile(inputFile,mesh.vertices,nAttribs,nMarkers)) return 1;
    }

    // Shift coordinates to near the origin to maximize floating-point precision.
    // Large coordinate offsets (e.g. x≈100000) cause orient2d/inCircle to lose
    // significant digits; centering on the bbox midpoint eliminates this.
    double shiftX=0, shiftY=0;
    if(!mesh.vertices.empty()){
        double minX=mesh.vertices[0].x, maxX=minX, minY=mesh.vertices[0].y, maxY=minY;
        for(auto& p:mesh.vertices){ minX=std::min(minX,p.x); maxX=std::max(maxX,p.x); minY=std::min(minY,p.y); maxY=std::max(maxY,p.y); }
        shiftX=(minX+maxX)*0.5; shiftY=(minY+maxY)*0.5;
        for(auto& p:mesh.vertices){ p.x-=shiftX; p.y-=shiftY; }
        for(auto& h:mesh.holes){ h.x-=shiftX; h.y-=shiftY; }
        for(auto& r:mesh.regions){ r.x-=shiftX; r.y-=shiftY; }
    }

    if(!opts.quiet){
        std::cerr<<"Input: "<<mesh.vertices.size()<<" vertices";
        if(!mesh.segments.empty()) std::cerr<<", "<<mesh.segments.size()<<" segments";
        if(!mesh.holes.empty()) std::cerr<<", "<<mesh.holes.size()<<" holes";
        if(!mesh.regions.empty()) std::cerr<<", "<<mesh.regions.size()<<" regions";
        std::cerr<<"\n";
    }

    auto tStart = std::chrono::high_resolution_clock::now();
    auto tPrev = tStart;
    auto elapsed = [&]() {
        auto now = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(now - tPrev).count();
        tPrev = now;
        return ms;
    };

    // 0. Optional PSLG cleanup
    if(opts.clean_pslg && opts.pslg){
        cleanPSLG(mesh, opts.clean_tol, opts.quiet);
        if(!opts.quiet) std::cerr<<"  PSLG cleanup: "<<elapsed()<<" ms\n";
    }

    // 1. Delaunay triangulation
    buildDelaunay(mesh);
    if(!opts.quiet) std::cerr<<"Initial Delaunay: "<<mesh.triangles.size()<<" triangles ("<<elapsed()<<" ms)\n";

    // 2. Enforce PSLG segments (CDT)
    mesh.rebuildAdjacency();
    if(opts.pslg){
        int miss = enforceConstraints(mesh, opts.quiet);

        // Split long unenforced segments at their midpoint and rebuild.
        // Only splits segments longer than a threshold; short segments that
        // fail are a robustness issue that splitting won't fix.
        for(int splitRound=0; splitRound<10 && miss>0; splitRound++){
            // Compute average segment length as a split threshold.
            // Only split segments significantly longer than the average
            // constrained segment, since splitting short failing segments
            // indicates a robustness issue rather than a length issue.
            double sumSegLen=0;
            for(auto& s:mesh.segments){
                auto& a=mesh.vertices[s.v0]; auto& b=mesh.vertices[s.v1];
                sumSegLen+=std::sqrt((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y));
            }
            double avgSegLen = sumSegLen / std::max(1,(int)mesh.segments.size());
            double minSplitLen = avgSegLen * 2.0;

            std::vector<std::vector<int>> v2t(mesh.vertices.size());
            for(int i=0;i<(int)mesh.triangles.size();i++){
                auto& t=mesh.triangles[i]; if(t.v[0]<0) continue;
                for(int j=0;j<3;j++) v2t[t.v[j]].push_back(i);
            }
            std::vector<int> toSplit;
            for(int i=0;i<(int)mesh.segments.size();i++){
                if(mesh.hasEdge(mesh.segments[i].v0,mesh.segments[i].v1,v2t)) continue;
                auto& a=mesh.vertices[mesh.segments[i].v0];
                auto& b=mesh.vertices[mesh.segments[i].v1];
                if(std::sqrt((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)) >= minSplitLen)
                    toSplit.push_back(i);
            }
            if(toSplit.empty()) break;

            int nSplit=0;
            for(int ii=(int)toSplit.size()-1; ii>=0; ii--){
                int si=toSplit[ii];
                auto& seg=mesh.segments[si];
                double mx=(mesh.vertices[seg.v0].x+mesh.vertices[seg.v1].x)/2;
                double my=(mesh.vertices[seg.v0].y+mesh.vertices[seg.v1].y)/2;
                int bm = std::max(mesh.vertices[seg.v0].marker,
                                  mesh.vertices[seg.v1].marker);
                int midIdx=(int)mesh.vertices.size();
                mesh.vertices.push_back({mx, my, midIdx, bm, {}});

                int v0=seg.v0, v1=seg.v1, mk=seg.marker;
                mesh.segments[si] = {v0, midIdx, mk};
                mesh.segments.insert(mesh.segments.begin()+si+1, {midIdx, v1, mk});
                nSplit++;
            }
            if(!opts.quiet)
                std::cerr<<"  Split "<<nSplit<<" long unenforced segments, rebuilding...\n";

            mesh.triangles.clear();
            buildDelaunay(mesh);
            mesh.rebuildAdjacency();
            miss = enforceConstraints(mesh, opts.quiet);
        }
    }
    if(!opts.quiet) std::cerr<<"  CDT: "<<elapsed()<<" ms\n";

    // 3. Remove holes via flood fill
    if(opts.pslg) removeHoles(mesh, opts);
    if(!opts.quiet) std::cerr<<"  Holes: "<<elapsed()<<" ms\n";

    // 4. Assign regions via flood fill on the coarse mesh
    if(opts.regions) assignRegions(mesh, opts);
    if(!opts.quiet) std::cerr<<"  Regions: "<<elapsed()<<" ms\n";

    // Note: region-0 triangles (gaps between domains) are kept in the mesh
    // during refinement to preserve topology. They are removed after refinement.

    int nInputVerts = (int)mesh.vertices.size();

    // 5. Quality refinement (new triangles inherit regions from containing triangle)
    if(opts.quality||opts.area_limit) refineQuality(mesh, opts, nInputVerts);
    if(!opts.quiet) std::cerr<<"  Refinement: "<<elapsed()<<" ms\n";

    // 7. Final cleanup: compact dead triangles, rebuild adjacency, extract edges
    {
        std::vector<Triangle> live;
        for(auto& t:mesh.triangles) if(t.v[0]>=0) live.push_back(t);
        mesh.triangles=live;
    }
    mesh.rebuildAdjacency();
    extractEdges(mesh);
    if(!opts.quiet) std::cerr<<"  Cleanup+edges: "<<elapsed()<<" ms\n";

    // ---- Convert PBC twin map to output pairs ----
    // The twin map was populated after reading .fem file (initial PBC node pairs)
    // and extended by synchronized segment splits during quality refinement.
    // Each entry appears twice (a→b and b→a); output each pair once.
    {
        std::set<std::pair<int,int>> seen;
        for(auto& [a, b] : mesh.pbc_twin){
            if(a==b) continue;
            auto key=std::make_pair(std::min(a,b), std::max(a,b));
            if(seen.insert(key).second){
                int type=mesh.pbc_node_type.count(a)?mesh.pbc_node_type[a]:0;
                mesh.pbc_pairs.push_back({a, b, type});
            }
        }
    }

    if(opts.jettison){
        std::vector<bool> used(mesh.vertices.size(),false);
        for(auto& t:mesh.triangles){used[t.v[0]]=used[t.v[1]]=used[t.v[2]]=true;}
        std::vector<int> vremap(mesh.vertices.size(),-1);
        std::vector<Point> newPts;
        for(int i=0;i<(int)mesh.vertices.size();i++){
            if(used[i]){vremap[i]=(int)newPts.size();newPts.push_back(mesh.vertices[i]);}
        }
        mesh.vertices=newPts;
        for(auto& t:mesh.triangles){t.v[0]=vremap[t.v[0]];t.v[1]=vremap[t.v[1]];t.v[2]=vremap[t.v[2]];}
        for(auto& e:mesh.edges){e.first=vremap[e.first];e.second=vremap[e.second];}
        for(auto& s:mesh.segments){s.v0=vremap[s.v0];s.v1=vremap[s.v1];}
        for(auto& p:mesh.pbc_pairs){p.node_a=vremap[p.node_a];p.node_b=vremap[p.node_b];}
    }

    if(!opts.quiet)
        std::cerr<<"Output: "<<mesh.vertices.size()<<" vertices, "<<mesh.triangles.size()<<" triangles\n";

    // Shift coordinates back to original position for output
    for(auto& p:mesh.vertices){ p.x+=shiftX; p.y+=shiftY; }

    std::string outBase=base+(opts.suppress_iter?"":".1");
    writeNodeFile(outBase+".node",mesh.vertices,nAttribs,nMarkers>0||opts.pslg,opts);
    writeEleFile(outBase+".ele",mesh,opts.regions,opts);
    if(opts.edges) writeEdgeFile(outBase+".edge",mesh,opts);
    if(opts.neighbors) writeNeighFile(outBase+".neigh",mesh,opts);
    if(opts.pslg&&!opts.no_poly_out) writePolyFile(outBase+".poly",mesh,opts);
    if(!mesh.pbc_pairs.empty()) writePbcFile(outBase+".pbc",mesh,opts);

    if(!opts.quiet){
        std::cerr<<"Wrote "<<outBase<<".node, "<<outBase<<".ele";
        if(opts.edges) std::cerr<<", "<<outBase<<".edge";
        if(opts.neighbors) std::cerr<<", "<<outBase<<".neigh";
        if(opts.pslg&&!opts.no_poly_out) std::cerr<<", "<<outBase<<".poly";
        if(!mesh.pbc_pairs.empty()) std::cerr<<", "<<outBase<<".pbc";
        std::cerr<<"\n";
    }
    return 0;
}
