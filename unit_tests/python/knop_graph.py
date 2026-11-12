#!/usr/bin/env python3
"""Layer-A transport dependency graph for the KNOP refresh (P1-T1h).

Builds the directed graph whose vertices are the active writable
radiation-store elements

    (channel, kc, ray, slot)     channel in {M,P}

where slot s corresponds to frequency s+1 and only slots
1..NCF-1 are writable (slot NCF is the never-written phantom).

Edges are derived from the source-level dependency rules traced in
docs/audits/KNOP_PRODUCTION_DEPENDENCY_GRAPH.md and re-verified
against plane/fg_j_cmf_nonmono_v1.f:

  * WRITE_STORE writes STORE(KC,LS,LS2-1)=FLD(JP2,LS2) with
    JP2=J_PNT(KS,LS), KC=COMP_DEPTH_IDX(KS): output (KC,LS,s) is the
    field at ray-point JP2, frequency s+1.
  * Compressing ray-point I (GAM<0) reads the lagged store at
    (KC1|KC2, LS, FREQ_IDX_ITER+NM) with KC1/KC2=COMP_DEPTH_IDX of the
    two depths bracketing RAY_PNT(I,LS)-ND_ADD, interpolated by
    W_INEXT on R_RAY vs R(K).  NM=NM_RAY(I) (1 when
    FG_FREQ_MATCH=F).
  * Expanding ray-point I (GAM>0) reads I_*_PREV = the previous
    frequency's field this sweep, giving a backward frequency reach.
  * The channel recurrences give the upwind source stencil:
    M marches I=1->NI (I_M(J) <- SOURCE_PRIME(1..J+1));
    P marches I=NI->1 (I_P(J) <- SOURCE_PRIME(J-1..NI)).

The expanding branch is handled by propagating field-level store
dependencies over frequency (deps_field[I][l]); the result is a
conservative structural graph that includes every dependency the
rules permit.  Finite-difference validation (a separate probe)
prunes it to the numerically-real edges.

This module is diagnostic-only: it reads KNOP_TOPO and writes graph
summaries; it never touches the solver.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from knop_topo import read_topo  # noqa: E402


# ----------------------------------------------------------------------
# Tarjan SCC (iterative) and Kahn topological sort.
# ----------------------------------------------------------------------
def tarjan_scc(n, adj):
    """adj: list of iterables of successor vertex ids (0-based).
    Returns list of SCCs, each a list of vertex ids."""
    index = [None] * n
    low = [0] * n
    onstk = [False] * n
    stk = []
    sccs = []
    counter = [0]

    for root in range(n):
        if index[root] is not None:
            continue
        work = [(root, iter(adj[root]))]
        index[root] = low[root] = counter[0]
        counter[0] += 1
        stk.append(root)
        onstk[root] = True
        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if index[w] is None:
                    index[w] = low[w] = counter[0]
                    counter[0] += 1
                    stk.append(w)
                    onstk[w] = True
                    work.append((w, iter(adj[w])))
                    advanced = True
                    break
                elif onstk[w]:
                    low[v] = min(low[v], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                u = work[-1][0]
                low[u] = min(low[u], low[v])
            if low[v] == index[v]:
                comp = []
                while True:
                    w = stk.pop()
                    onstk[w] = False
                    comp.append(w)
                    if w == v:
                        break
                sccs.append(comp)
    return sccs


def topo_order(n, adj):
    """Kahn topological order of a DAG; returns order or None if cyclic."""
    indeg = [0] * n
    radj = [[] for _ in range(n)]
    for v in range(n):
        for w in adj[v]:
            indeg[w] += 1
            radj[w].append(v)
    from collections import deque
    q = deque([v for v in range(n) if indeg[v] == 0])
    order = []
    while q:
        v = q.popleft()
        order.append(v)
        for w in adj[v]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    return order if len(order) == n else None


# ----------------------------------------------------------------------
# Graph construction.
# ----------------------------------------------------------------------
def build_layerA(topo, nm_offset=None):
    """Build the conservative Layer-A transport graph.

    Returns (vertices, edges, meta) where vertices is a list of
    (channel,kc,ray,slot) tuples (1-based slot/freq semantics kept in
    meta), edges is a list of (src_vertex_id, dst_vertex_id), and
    meta carries per-ray/channel diagnostics.

    Vertices use channel index 0=M,1=P.  Only compressing depths with
    a valid J_PNT ray-point on the ray produce output vertices.
    """
    a = topo["arrays"]
    dims = topo["dims"]
    ND, NP, NCF, NCOMP = (dims["KN_ND"], dims["KN_NP"], dims["KN_NCF"],
                          dims["N_COMP"])
    NRMAX = dims["NRAY_MAX"]
    ND_ADD = topo["icnf"]["ND_ADD"]

    COMP_DEPTH_IDX = a["COMP_DEPTH_IDX"]          # (ND,) -> kc or 0
    J_PNT = a["J_PNT"]                            # (ND,NP)
    RAY_PNT = a["RAY_PNT"]                        # (NRMAX,NP)
    NI_RAY = a["NI_RAY"]                          # (NP,)
    GAM = a["GAM"]                                # (NRMAX,NP)
    KN_FREQ = a["KN_FREQ"]                        # (NCF,)

    fg_match = topo["lcfg"].get("FG_FREQ_MATCH", False)
    if nm_offset is None:
        # FG_FREQ_MATCH=F -> NM_RAY=1.  For =T the per-(ray,point,freq)
        # offset is a derived quantity; use the conservative max reach
        # supported by the frequency grid unless a caller overrides.
        nm_offset = 1 if not fg_match else max(1, NCF // 2)

    # compressing store depths: KS (1..ND) -> kc=COMP_DEPTH_IDX(KS)>=1
    comp_KS = [KS for KS in range(ND) if COMP_DEPTH_IDX[KS] >= 1]

    # slot range: writable store slots 1..NCF-1 (slot NCF phantom)
    S_LO, S_HI = 1, NCF - 1

    # vertices[(c,kc,LS,s)] -> id ; id -> tuple
    vid = {}
    verts = []

    def vkey(c, kc, ls, s):
        return (c, int(kc), int(ls), int(s))

    def get_vid(c, kc, ls, s):
        k = vkey(c, kc, ls, s)
        if k not in vid:
            vid[k] = len(verts)
            verts.append(k)
        return vid[k]

    # Pre-create all output vertices so edge endpoints always exist.
    for c in (0, 1):
        for ls in range(1, NP + 1):
            ni = int(NI_RAY[ls - 1])
            for KS in comp_KS:
                kc = int(COMP_DEPTH_IDX[KS])
                jp2 = int(J_PNT[KS - 1, ls - 1])
                if not (1 <= jp2 <= ni):
                    continue
                for s in range(S_LO, S_HI + 1):
                    get_vid(c, kc, ls, s)

    edges = set()

    # Per (channel, ray): field-level store-dependency propagation.
    # deps_field[I][l] : set of (kc_in, slot_in) read by field at
    # ray-point I (1-based) during frequency l (1-based, l=1..NCF).
    # slot_in is an input store slot; only 1..NCF-1 are real vertices.
    for c, march in ((0, "M"), (1, "P")):
        for ls in range(1, NP + 1):
            ni = int(NI_RAY[ls - 1])
            if ni < 1:
                continue

            def upwind(i):
                # 1-based ray-point stencil (documented conservative).
                if march == "M":
                    return range(1, min(i + 1, ni) + 1)
                return range(max(1, i - 1), ni + 1)

            # deps_field[I][l], I=1..ni, l=1..NCF
            deps_field = [[set() for _ in range(NCF + 1)]
                          for _ in range(ni + 1)]

            for l in range(1, NCF + 1):
                for i in range(1, ni + 1):
                    acc = set()
                    for ip in upwind(i):
                        g = GAM[ip - 1, ls - 1]
                        if g < 0.0:
                            # compressing donor read at slot l
                            k = int(RAY_PNT[ip - 1, ls - 1]) - ND_ADD
                            donors = set()
                            for kd in (k, k + 1):
                                if 1 <= kd <= ND:
                                    kc_in = int(COMP_DEPTH_IDX[kd - 1])
                                    if kc_in >= 1:
                                        donors.add(kc_in)
                            slot_in = l  # FREQ_IDX_ITER(l-1)+NM(1)=l
                            if S_LO <= slot_in <= S_HI:
                                for kc_in in donors:
                                    acc.add((kc_in, slot_in))
                        elif g > 0.0:
                            if l - 1 >= 1:
                                acc |= deps_field[ip][l - 1]
                        # g==0: no coupling
                    deps_field[i][l] = acc

            # project field deps onto store-output vertices
            for KS in comp_KS:
                kc = int(COMP_DEPTH_IDX[KS])
                jp2 = int(J_PNT[KS - 1, ls - 1])
                if not (1 <= jp2 <= ni):
                    continue
                for s in range(S_LO, S_HI + 1):
                    l = s + 1  # field frequency that writes slot s
                    dst = vid[vkey(c, kc, ls, s)]
                    for (kc_in, slot_in) in deps_field[jp2][l]:
                        src = vid.get(vkey(c, kc_in, ls, slot_in))
                        if src is not None:
                            edges.add((src, dst))

    edges = sorted(edges)
    adj = [[] for _ in range(len(verts))]
    radj = [[] for _ in range(len(verts))]
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)

    meta = {
        "ncomp": NCOMP, "ncf": NCF, "nd": ND, "np": NP,
        "comp_KS": comp_KS, "nm_offset": nm_offset,
        "fg_freq_match": fg_match,
    }
    return verts, edges, adj, radj, meta


def graph_stats(verts, edges, adj):
    n = len(verts)
    self_loops = sum(1 for u, v in edges if u == v)
    sccs = tarjan_scc(n, adj)
    nontriv = [c for c in sccs if len(c) > 1]
    self_scc = [c for c in sccs if len(c) == 1 and
                any(w == c[0] for w in adj[c[0]])]
    order = topo_order(n, adj)
    # longest chain (only meaningful if DAG)
    maxchain = None
    if order is not None:
        depth = [0] * n
        for v in order:
            for w in adj[v]:
                depth[w] = max(depth[w], depth[v] + 1)
        maxchain = max(depth) if depth else 0
    return {
        "n_vertices": n,
        "n_edges": len(edges),
        "self_loops": self_loops,
        "n_scc": len(sccs),
        "n_nontrivial_scc": len(nontriv) + len(self_scc),
        "largest_scc": max((len(c) for c in sccs), default=0),
        "is_dag": order is not None,
        "topo_order": order,
        "max_chain": maxchain,
        "sccs": sccs,
    }


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "KNOP_TOPO"
    topo = read_topo(path)
    verts, edges, adj, radj, meta = build_layerA(topo)
    st = graph_stats(verts, edges, adj)
    print(f"vertices={st['n_vertices']} edges={st['n_edges']} "
          f"selfloops={st['self_loops']}")
    print(f"scc={st['n_scc']} nontrivial={st['n_nontrivial_scc']} "
          f"largest={st['largest_scc']} dag={st['is_dag']} "
          f"maxchain={st['max_chain']}")
    # channel/depth/slot histogram of edges
    from collections import Counter
    chan = Counter()
    slotdiff = Counter()
    for u, v in edges:
        cu, _, _, su = verts[u]
        cv, _, _, sv = verts[v]
        chan[(cu, cv)] += 1
        slotdiff[sv - su] += 1
    print("channel pairs:", dict(chan))
    print("slot diffs (dst-src):", dict(sorted(slotdiff.items())))


if __name__ == "__main__":
    main()
