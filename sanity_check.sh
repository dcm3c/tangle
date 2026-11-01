#!/bin/bash
# sanity_check.sh -- run tangle over the full test battery and flag any mesh that
# is INVALID (inverted or degenerate triangles) or fails to terminate. Run this
# after any mesh-affecting change to confirm nothing regressed.
#
#   ./sanity_check.sh                 # ./tangle.exe (or ./tangle), -pq30, 60s timeout
#   ANGLE=33 ./sanity_check.sh        # different quality angle
#   TANGLE=./other.exe ./sanity_check.sh
#
# This checks VALIDITY + TERMINATION only (the "did I break it" guard). It does
# NOT check byte-identity against a reference -- that's a separate, stricter test.
# Known residual NOT caught here: test_orig at -pq33aA (FEMM path) still caps the
# Steiner budget; see CLAUDE.md.

cd "$(dirname "$0")" || exit 1

TANGLE="${TANGLE:-./tangle.exe}"
[ -x "$TANGLE" ] || TANGLE=./tangle
TIMEOUT="${TIMEOUT:-60}"
ANGLE="${ANGLE:-30}"

# .poly battery: normal + hard + ex/ + test_orig. Missing files are skipped.
# (ex/box is deliberately omitted: its 'outer box' is a single segment, so it
#  encloses no area -- Triangle and tangle both produce an empty mesh, and
#  tangle now reports that as an error rather than writing empty output.)
POLY="myExample rmb motor1 motor2 roters1b CircleOnBoundary bubbles dslim1 \
bigGrid2 grid2 grid3 grid4 lrk two_circles wedge_pbc adjpbc nearline vf test_orig \
ex/A ex/airfoil_exterior ex/airfoil_interior ex/double_hex ex/double_hex3 ex/la"

# .fem battery: native FEMM reader; outputs <base>.node/.ele (no .1.).
# 2horseAGE is the PBC+AGE exerciser (5 anti-periodic pairs + air gap element).
FEM="Antunes lrk my_SPM_motor roters1b bigModel TinyTriAngleTest 2horseAGE 2horseAGE_dup kelvin_half kelvin_half_auxfirst"

# Other FEMM physics (heat .feh, electrostatics .fee): BdryFormat numbering for
# periodic/anti-periodic differs per physics type and nothing else covers it.
OTHER="heat_rad.feh heat_per.feh heat_anti.feh strips.fee"

# Inputs that MUST be refused, as input:expected-exit pairs. pbc_mismatch: a
# periodic boundary whose two sides cannot be paired is a mesh error (4), never a
# silent natural boundary. badseg_idx / badarc_idx: a segment or arc naming a
# vertex outside the vertex list is a parse error (3); before the check it was
# undefined behavior (segfault or runaway memory in the arc chordizer).
NEG_POLY="pbc_mismatch:4 badseg_idx:3 badarc_idx:3"

# Fast O(n) validity check. Prints "<tris> <inverted> <degenerate>" or "ERR".
chk() {  # $1=node $2=ele
  python - "$1" "$2" <<'PY'
import sys
try:
    V={}
    L=open(sys.argv[1]).read().split('\n'); nn=int(L[0].split()[0]); i=1; c=0
    while c<nn:
        t=L[i].split(); i+=1
        if not t or t[0].startswith('#'): continue
        V[int(t[0])]=(float(t[1]),float(t[2])); c+=1
    L=open(sys.argv[2]).read().split('\n'); ne=int(L[0].split()[0]); i=1; c=0
    inv=deg=0
    while c<ne:
        t=L[i].split(); i+=1
        if not t or t[0].startswith('#'): continue
        a,b,cc=int(t[1]),int(t[2]),int(t[3]); c+=1
        p,q,r=V[a],V[b],V[cc]
        ar=(q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])
        if ar<-1e-12: inv+=1
        elif abs(ar)<1e-12: deg+=1
    print(f"{ne} {inv} {deg}")
except Exception:
    print("ERR")
PY
}

pass=0; fail=0; skip=0

run_case() {  # $1=input-arg  $2=flags  $3=node  $4=ele  $5=label
  local arg="$1" flags="$2" node="$3" ele="$4" label="$5"
  rm -f "$node" "$ele" "${node%.node}.pbc"   # stale .pbc would let a pairing check pass vacuously
  local t0 t1 e
  t0=$(date +%s.%N)
  timeout "$TIMEOUT" $TANGLE $flags "$arg" >/dev/null 2>&1; e=$?
  t1=$(date +%s.%N)
  local dt; dt=$(awk "BEGIN{printf \"%.2f\", $t1-$t0}")
  if [ $e -eq 124 ]; then
    printf "%-26s %10s %6s %6s  %6ss  *** HANG ***\n" "$label" - - - "$dt"; fail=$((fail+1)); return
  fi
  if [ $e -ne 0 ]; then
    printf "%-26s %10s %6s %6s  %6ss  *** exit=%d ***\n" "$label" - - - "$dt" "$e"; fail=$((fail+1)); return
  fi
  if [ ! -f "$ele" ]; then
    printf "%-26s %10s %6s %6s  %6ss  *** NO OUTPUT ***\n" "$label" - - - "$dt"; fail=$((fail+1)); return
  fi
  local out; out=$(chk "$node" "$ele")
  local tris inv deg; read tris inv deg <<< "$out"
  if [ "$tris" = "ERR" ]; then
    printf "%-26s %10s %6s %6s  %6ss  *** check error ***\n" "$label" - - - "$dt"; fail=$((fail+1)); return
  fi
  if [ "${inv:-0}" != "0" ] || [ "${deg:-0}" != "0" ]; then
    printf "%-26s %10s %6s %6s  %6ss  <<< FAIL (invalid)\n" "$label" "$tris" "$inv" "$deg" "$dt"; fail=$((fail+1))
  else
    printf "%-26s %10s %6s %6s  %6ss  ok\n" "$label" "$tris" "$inv" "$deg" "$dt"; pass=$((pass+1))
  fi
}

echo "tangle = $TANGLE   |   -pq${ANGLE}   |   timeout ${TIMEOUT}s"
printf "%-26s %10s %6s %6s  %7s  %s\n" "case" "tris" "inv" "deg" "time" "status"
printf -- "----------------------------------------------------------------------\n"

for f in $POLY; do
  [ -f "${f}.poly" ] || { skip=$((skip+1)); continue; }
  run_case "$f" "-pq${ANGLE} -e -Q" "${f}.1.node" "${f}.1.ele" "$f"
done
for f in $FEM; do
  [ -f "${f}.fem" ] || { skip=$((skip+1)); continue; }
  run_case "${f}.fem" "" "${f}.node" "${f}.ele" "${f}.fem"
done
for f in $OTHER; do
  [ -f "$f" ] || { skip=$((skip+1)); continue; }
  run_case "$f" "" "${f%.*}.node" "${f%.*}.ele" "$f"
done
for fx in $NEG_POLY; do
  f=${fx%%:*}; want=${fx##*:}
  [ -f "${f}.poly" ] || { skip=$((skip+1)); continue; }
  rm -f "${f}.1.node" "${f}.1.ele"
  timeout "$TIMEOUT" $TANGLE -pq${ANGLE} -Q "${f}.poly" >/dev/null 2>&1; e=$?
  if [ $e -eq "$want" ] && [ ! -f "${f}.1.node" ]; then
    printf "%-26s %10s %6s %6s  %7s  ok (refused, exit %d)\n" "$f" - - - - "$e"; pass=$((pass+1))
  else
    extra=""; [ -f "${f}.1.node" ] && extra=", output written"
    printf "%-26s %10s %6s %6s  %7s  <<< FAIL (expected exit %d, got %d%s)\n" "$f" - - - - "$want" "$e" "$extra"; fail=$((fail+1))
  fi
done

# ---- Pairing / periodic-boundary checks on the outputs produced above ----
# These guard bugs that leave a VALID mesh behind (so the checks above pass)
# while the .pbc is wrong: the AGE ring shift after a node merge, and the
# Kelvin-corner desync that silently unpaired an anti-periodic diameter.
check() {  # $1=label $2=0/1 ok $3=detail
  if [ "$2" = 1 ]; then printf "%-26s %s  ok
" "$1" "$3"; pass=$((pass+1))
  else printf "%-26s %s  <<< FAIL
" "$1" "$3"; fail=$((fail+1)); fi
}
# .pbc summary: "<pairs> <types> <axis-nodes> <axis-paired>" (axis = nodes on y=0)
pbcinfo() {  # $1=node $2=pbc (may be missing)
  python - "$1" "$2" <<'PY2'
import sys,os
L=open(sys.argv[1]).read().splitlines(); n=int(L[0].split()[0])
ax={int(L[i].split()[0]) for i in range(1,n+1) if abs(float(L[i].split()[2]))<1e-9}
pr=set(); ty=set(); m=0
if os.path.exists(sys.argv[2]):
    P=open(sys.argv[2]).read().splitlines(); m=int(P[0])
    for i in range(1,m+1):
        t=P[i].split(); pr.add(int(t[1])); pr.add(int(t[2])); ty.add(t[3])
print(m, "+".join(sorted(ty)) or "-", len(ax), len(ax&pr))
PY2
}
printf -- "----------------------------------------------------------------------
"
if [ -f 2horseAGE.pbc ] && [ -f 2horseAGE_dup.pbc ]; then
  ok=1; for e in node ele pbc; do cmp -s "2horseAGE.$e" "2horseAGE_dup.$e" || ok=0; done
  check "2horseAGE_dup == 2horseAGE" $ok "node/ele/pbc byte-identical after merging the duplicate point"
fi
for f in kelvin_half kelvin_half_auxfirst; do
  [ -f "$f.node" ] || continue
  read m ty ax axp <<< "$(pbcinfo $f.node $f.pbc)"
  [ "$ax" -gt 0 ] && [ "$ax" = "$axp" ] && ok=1 || ok=0
  check "$f axis pairing" $ok "$axp of $ax y=0 nodes paired, $m pairs ($ty)"
done
if [ -f kelvin_half.pbc ] && [ -f kelvin_half_auxfirst.pbc ]; then
  ok=1; cmp -s kelvin_half.ele kelvin_half_auxfirst.ele || ok=0; cmp -s kelvin_half.pbc kelvin_half_auxfirst.pbc || ok=0
  check "kelvin_half order-indep." $ok "same .ele and .pbc whichever boundary property comes first"
fi
heatcheck() {  # $1=base $2=expected types ("-" = no .pbc at all)
  [ -f "$1.node" ] || return
  read m ty ax axp <<< "$(pbcinfo $1.node $1.pbc)"
  if [ "$2" = "-" ]; then [ "$m" = 0 ] && ok=1 || ok=0; else [ "$m" -gt 0 ] && [ "$ty" = "$2" ] && ok=1 || ok=0; fi
  check "$1 pairing" $ok "$m pairs, types $ty (expected $2)"
}
heatcheck heat_rad  -    # radiation (BdryFormat 3) is NOT periodic
heatcheck heat_per  0    # heat periodic is 4
heatcheck heat_anti 1    # heat anti-periodic is 5
heatcheck strips    0    # .fee periodic is 3

# .pbc node numbers must index the .node file written alongside it (same base,
# -z or not). wedge_pbc is a periodic ROTATION about the origin, so every pair
# must join two nodes of equal radius; an off-by-one in the .pbc numbering (it
# was always 0-based while .node is 1-based without -z) breaks that at once.
if [ -f wedge_pbc.1.node ] && [ -f wedge_pbc.1.pbc ]; then
  res=$(python - wedge_pbc.1.node wedge_pbc.1.pbc <<'PY3'
import sys, math
L = open(sys.argv[1]).read().splitlines(); n = int(L[0].split()[0])
V = {}
for r in L[1:n+1]:
    t = r.split(); V[int(t[0])] = (float(t[1]), float(t[2]))
P = open(sys.argv[2]).read().splitlines(); m = int(P[0]); bad = 0
for r in P[1:m+1]:
    t = r.split(); a, b = int(t[1]), int(t[2])
    if a not in V or b not in V or abs(math.hypot(*V[a]) - math.hypot(*V[b])) > 1e-6 * math.hypot(*V[a]):
        bad += 1
print(m, bad)
PY3
)
  read m bad <<< "$res"
  [ "${bad:-1}" = 0 ] && [ "${m:-0}" -gt 0 ] && ok=1 || ok=0
  check "wedge_pbc .pbc numbering" $ok "$m pairs, $bad not equal-radius in .node numbering"
fi

printf -- "----------------------------------------------------------------------\n"
printf "PASS=%d  FAIL=%d  SKIP=%d\n" "$pass" "$fail" "$skip"
[ "$fail" -eq 0 ]
