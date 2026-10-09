# -*- coding: utf-8 -*-
"""通用拉伸体 STL 重建器 v2 —— 轴向自适应 + 多内孔

从任意「平板拉伸件」STL 自动识别拉伸轴、提取截面轮廓 → 重新三角化 → 输出绝对干净的网格。
消除 CAD 导出残留的非法边、未缝合顶点、退化面。

v2 改进（相对 v1）：
  * 自动检测拉伸轴（取尺寸最小的轴），不再假设厚度一定沿 Z
  * 支持任意数量内孔：四方向射线 + 桥接线可见性检查，一次性构造带全部切口的单环
  * 输出时保持原坐标系方位

用法: python rebuild_extrude.py <输入.stl> <输出.stl> [厚度mm]
      厚度为 0 或省略时沿用原文件厚度
"""
import struct, sys, math
from collections import defaultdict

Q = 1e-4
EPS = 1e-9

def load(path):
    data = open(path, 'rb').read()
    n = struct.unpack('<I', data[80:84])[0]
    n2 = (len(data) - 84) // 50
    n = n if n == n2 else n2
    tris, off = [], 84
    for i in range(n):
        v = struct.unpack('<12fH', data[off:off+50])
        tris.append((v[3:6], v[6:9], v[9:12])); off += 50
    return tris

def detect_axis(tris):
    """拉伸轴 = 尺寸最小的轴；返回 (轴索引, 各轴最小值, 各轴最大值)"""
    rng, lo, hi = [], [], []
    for k in range(3):
        vs = [p[k] for t in tris for p in t]
        rng.append(max(vs) - min(vs)); lo.append(min(vs)); hi.append(max(vs))
    return rng.index(min(rng)), lo, hi

def extract_rings(tris, w, wmin):
    """提取拉伸轴起点平面上的 2D 轮廓环（闭合环，首尾重复）"""
    uv = [i for i in range(3) if i != w]
    def flat(p): return (p[uv[0]], p[uv[1]])
    plane = [t for t in tris if all(abs(p[w] - wmin) < 1e-4 for p in t)]
    def vk2(p): return (round(p[0]/Q), round(p[1]/Q))
    dirs = defaultdict(list)
    for t in plane:
        a, b, c = vk2(flat(t[0])), vk2(flat(t[1])), vk2(flat(t[2]))
        for u, v in ((a, b), (b, c), (c, a)):
            if u != v: dirs[(u, v)].append(t)
    boundary = [(u, v) for (u, v) in dirs if (v, u) not in dirs]
    nxt = defaultdict(list)
    for u, v in boundary: nxt[u].append(v)
    rings, visited = [], set()
    for u, v in boundary:
        if (u, v) in visited: continue
        ring, cur, guard = [u], u, 0
        while guard < 20000:
            guard += 1
            if not nxt[cur] or (cur, nxt[cur][0]) in visited: break
            nv = nxt[cur][0]; visited.add((cur, nv))
            if nv == ring[0]: ring.append(nv); break
            if nv in ring: break
            ring.append(nv); cur = nv
        if len(ring) > 3: rings.append([(p[0]*Q, p[1]*Q) for p in ring])
    return rings

def signed_area(p):
    s = 0.0
    for i in range(len(p)-1):
        s += p[i][0]*p[i+1][1] - p[i+1][0]*p[i][1]
    return s / 2.0

def point_in_tri(p, a, b, c):
    def sgn(a, b, c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    d1, d2, d3 = sgn(a, b, p), sgn(b, c, p), sgn(c, a, p)
    E = 1e-9
    if abs(d1) < E or abs(d2) < E or abs(d3) < E: return False
    neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (neg and pos)

def seg_seg_hit(p1, p2, p3, p4):
    """线段严格相交（端点接触 / 共线重叠 均不算）"""
    def s(a, b, c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    d1, d2 = s(p3, p4, p1), s(p3, p4, p2)
    d3, d4 = s(p1, p2, p3), s(p1, p2, p4)
    E = 1e-9
    return (((d1 > E and d2 < -E) or (d1 < -E and d2 > E)) and
            ((d3 > E and d4 < -E) or (d3 < -E and d4 > E)))

def bridge_clear(P, Qp, outer, oi, holes):
    """桥接线 P-Qp 不穿越外环其他边、也不穿越任何孔"""
    n = len(outer) - 1
    for j in range(n):
        if j == oi: continue
        if seg_seg_hit(P, Qp, outer[j], outer[j+1]): return False
    for h in holes:
        m = len(h) - 1
        for j in range(m):
            if seg_seg_hit(P, Qp, h[j], h[j+1]): return False
    return True

def find_bridge(outer, hole, other_holes):
    """为单个孔找最短且无穿越的桥接线。
    遍历孔的所有顶点 × 外环所有边 × 水平/竖直射线，返回 (长度, 外环边索引, 交点Q, 孔顶点索引)"""
    n, m = len(outer)-1, len(hole)-1
    best = None
    blockers = [hole] + list(other_holes)
    for k in range(m):
        P = hole[k]
        for i in range(n):
            a, b = outer[i], outer[i+1]
            # 水平射线（沿 ±u）
            if abs(a[1]-b[1]) >= EPS and min(a[1], b[1])-1e-9 <= P[1] <= max(a[1], b[1])+1e-9:
                t = (P[1]-a[1])/(b[1]-a[1]); x = a[0] + t*(b[0]-a[0])
                if abs(x - P[0]) > 1e-9:
                    Qp = (x, P[1]); L = abs(x - P[0])
                    if (best is None or L < best[0]) and bridge_clear(P, Qp, outer, i, blockers):
                        best = (L, i, Qp, k, t)
            # 竖直射线（沿 ±v）
            if abs(a[0]-b[0]) >= EPS and min(a[0], b[0])-1e-9 <= P[0] <= max(a[0], b[0])+1e-9:
                t = (P[0]-a[0])/(b[0]-a[0]); y = a[1] + t*(b[1]-a[1])
                if abs(y - P[1]) > 1e-9:
                    Qp = (P[0], y); L = abs(y - P[1])
                    if (best is None or L < best[0]) and bridge_clear(P, Qp, outer, i, blockers):
                        best = (L, i, Qp, k, t)
    return best

def triangulate(poly):
    """耳切法三角化（poly 为闭合环，首尾可不重复）"""
    pts = [tuple(p) for p in poly]
    if len(pts) > 1 and abs(pts[0][0]-pts[-1][0]) < 1e-12 and abs(pts[0][1]-pts[-1][1]) < 1e-12:
        pts = pts[:-1]
    idx = list(range(len(pts))); tris = []
    guard = 0
    while len(idx) > 3 and guard < 200000:
        guard += 1; cut = False
        for i in range(len(idx)):
            ia, ib, ic = idx[i-1], idx[i], idx[(i+1) % len(idx)]
            a, b, c = pts[ia], pts[ib], pts[ic]
            cr = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
            if cr <= 1e-7: continue
            ok = True
            for j in idx:
                if j in (ia, ib, ic): continue
                if point_in_tri(pts[j], a, b, c): ok = False; break
            if ok:
                tris.append((a, b, c)); idx.pop(i); cut = True; break
        if not cut: break
    if len(idx) == 3: tris.append((pts[idx[0]], pts[idx[1]], pts[idx[2]]))
    elif len(idx) > 3: print("   [警告] 三角化未完成，剩余 %d 个顶点" % len(idx))
    return tris

def normal_of(t):
    a, b, c = t
    u = (b[0]-a[0], b[1]-a[1], b[2]-a[2]); v = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
    n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
    L = math.sqrt(n[0]**2+n[1]**2+n[2]**2)
    return (n[0]/L, n[1]/L, n[2]/L) if L > 1e-12 else (0.0, 0.0, 1.0)

def close_ring(r):
    if len(r) > 1 and (abs(r[0][0]-r[-1][0]) > 1e-9 or abs(r[0][1]-r[-1][1]) > 1e-9):
        return list(r) + [r[0]]
    return list(r)

def build(rings, T, w, wmin):
    uv = [i for i in range(3) if i != w]
    def make3(u, v, wv):
        p = [0.0, 0.0, 0.0]
        p[uv[0]] = u; p[uv[1]] = v; p[w] = wv
        return (p[0], p[1], p[2])

    rings = [close_ring(r) for r in rings]
    areas = [abs(signed_area(r)) for r in rings]
    mi = areas.index(max(areas))
    outer = [tuple(p) for p in rings[mi]]
    holes = [[tuple(p) for p in r] for i, r in enumerate(rings) if i != mi]
    if signed_area(outer) < 0: outer = outer[::-1]
    for i, h in enumerate(holes):
        if signed_area(h) > 0: holes[i] = h[::-1]

    # ---- 为每个孔规划桥接线 ----
    qmap = defaultdict(list)   # 外环边索引 -> [(交点参数t, 交点Q, 孔环, 孔顶点索引)]
    for idx, h in enumerate(holes):
        others = [holes[j] for j in range(len(holes)) if j != idx]
        got = find_bridge(outer, h, others)
        if got is None:
            raise RuntimeError("内孔 #%d（%.2f mm²）未找到无穿越的桥接位置" % (idx+1, areas[idx]))
        L, oi, Qp, k, t = got
        qmap[oi].append((t, Qp, h, k))

    # ---- 一次性构造带全部切口的单环 ----
    n = len(outer) - 1
    ring = []
    for i in range(n):
        ring.append(outer[i])
        if i in qmap:
            for t, Qp, h, k in sorted(qmap[i], key=lambda z: z[0]):
                m = len(h) - 1
                P = h[k]
                ring.append(Qp); ring.append(P)
                for j in range(m): ring.append(h[(k+j) % m])
                ring.append(P); ring.append(Qp)
    ring.append(outer[n])
    ded = []
    for p in ring:
        if not ded or abs(p[0]-ded[-1][0]) > 1e-9 or abs(p[1]-ded[-1][1]) > 1e-9:
            ded.append(p)
    if len(ded) > 1 and abs(ded[0][0]-ded[-1][0]) < 1e-9 and abs(ded[0][1]-ded[-1][1]) < 1e-9:
        ded.pop()
    face2d = triangulate(ded)

    # ---- 侧壁环：外环（插入所有切口点）+ 各孔原始环 ----
    side_outer = []
    for i in range(n):
        side_outer.append(outer[i])
        if i in qmap:
            for t, Qp, h, k in sorted(qmap[i], key=lambda z: z[0]):
                side_outer.append(Qp)
    side_outer.append(outer[n])
    side_rings = [side_outer] + [list(h) for h in holes]

    tris = []
    for (a, b, c) in face2d:
        tris.append([make3(a[0], a[1], wmin+T), make3(b[0], b[1], wmin+T), make3(c[0], c[1], wmin+T)])
    for (a, b, c) in face2d:
        tris.append([make3(c[0], c[1], wmin), make3(b[0], b[1], wmin), make3(a[0], a[1], wmin)])
    for ring2 in side_rings:
        m = len(ring2)
        for i in range(m):
            a, b = ring2[i], ring2[(i+1) % m]
            if math.hypot(b[0]-a[0], b[1]-a[1]) < 1e-9: continue
            nrm = ((b[1]-a[1]), -(b[0]-a[0]))
            L = math.hypot(*nrm); nrm = (nrm[0]/L, nrm[1]/L)
            for t2 in ([make3(a[0], a[1], wmin), make3(b[0], b[1], wmin), make3(b[0], b[1], wmin+T)],
                       [make3(a[0], a[1], wmin), make3(b[0], b[1], wmin+T), make3(a[0], a[1], wmin+T)]):
                nn = normal_of(t2)
                if nn[uv[0]]*nrm[0] + nn[uv[1]]*nrm[1] < 0: t2 = [t2[0], t2[2], t2[1]]
                tris.append(t2)
    return tris, face2d, areas

def volume(tris):
    s = 0.0
    for a, b, c in tris:
        s += (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
              + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0
    return s

def check(tris):
    def vk(p): return (round(p[0]/Q), round(p[1]/Q), round(p[2]/Q))
    def ek(a, b):
        ka, kb = vk(a), vk(b); return (ka, kb) if ka <= kb else (kb, ka)
    e = defaultdict(int)
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            if vk(a) == vk(b): return -1, "存在退化三角形"
            e[ek(a, b)] += 1
    return sum(1 for v in e.values() if v != 2), None

def save_stl(path, tris):
    out = bytearray()
    hdr = b'Rebuilt by WorkBuddy STL Toolkit v2'
    out += hdr + b' ' * (80-len(hdr))
    out += struct.pack('<I', len(tris))
    for t in tris:
        out += struct.pack('<3f', *normal_of(t))
        for p in t: out += struct.pack('<3f', *p)
        out += struct.pack('<H', 0)
    open(path, 'wb').write(bytes(out))

def main():
    src, dst = sys.argv[1], sys.argv[2]
    T = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    tris_in = load(src)
    w, lo, hi = detect_axis(tris_in)
    print("拉伸轴: %s  (X=%.2f Y=%.2f Z=%.2f)" % ("XYZ"[w], hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]))
    if T <= 0: T = hi[w] - lo[w]
    print("原厚度 %.3f mm -> 目标厚度 %.3f mm" % (hi[w]-lo[w], T))
    rings = extract_rings(tris_in, w, lo[w])
    minu = min(p[0] for r in rings for p in r)
    minv = min(p[1] for r in rings for p in r)
    rings = [[(p[0]-minu, p[1]-minv) for p in r] for r in rings]
    areas = [abs(signed_area(r)) for r in rings]
    mi = areas.index(max(areas))
    print("环: 主体 %.2f mm2 + %d 个内孔 [%s] = 净 %.2f mm2"
          % (areas[mi], len(areas)-1,
             ", ".join("%.2f" % a for i, a in enumerate(areas) if i != mi),
             areas[mi] - (sum(areas) - areas[mi])))
    tris, f2d, areas = build(rings, T, w, lo[w])
    vol = volume(tris)
    bad, err = check(tris)
    a2d = sum(abs((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))/2 for a, b, c in f2d)
    print("2D 三角化 %d 个 | 净面积 %.2f mm2 | 三角形总数 %d | 体积 %.2f mm3"
          % (len(f2d), a2d, len(tris), abs(vol)))
    if err: print("!! %s" % err)
    else: print("非流形边: %d %s" % (bad, "(合格)" if bad == 0 else "(异常!)"))
    save_stl(dst, tris)
    print("已保存: %s" % dst)

if __name__ == '__main__':
    main()
