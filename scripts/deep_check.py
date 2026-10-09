# -*- coding: utf-8 -*-
"""STL 深度检查：共面重叠面 / 自相交 / 薄壁 / 短边 / 近重合顶点 / 壳体数"""
import struct, sys, math, itertools
from collections import defaultdict

path = sys.argv[1]
data = open(path, 'rb').read()
n = struct.unpack('<I', data[80:84])[0]
n2 = (len(data) - 84) // 50
n = n if n == n2 else n2
tris = []
off = 84
for i in range(n):
    v = struct.unpack('<12fH', data[off:off+50])
    tris.append((v[3:6], v[6:9], v[9:12]))
    off += 50
print("三角形数: %d" % len(tris))

def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a, b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def norm(a):
    L = math.sqrt(dot(a, a))
    return (a[0]/L, a[1]/L, a[2]/L) if L > 0 else (0, 0, 0)

tri_norm = []
for t in tris:
    tri_norm.append(norm(cross(sub(t[1], t[0]), sub(t[2], t[0]))))

# ---------- 1. 共面重叠面（同一平面 + 面积重叠）----------
def tri_area(t):
    return 0.5 * math.sqrt(dot(cross(sub(t[1], t[0]), sub(t[2], t[0])),
                               cross(sub(t[1], t[0]), sub(t[2], t[0]))))
print("\n=== 1. 共面且重叠的三角形对（真正的重叠面 = 非法）===")
coplanar_dup = 0
DOTS = 1 - 1e-6
for i, j in itertools.combinations(range(len(tris)), 2):
    ni, nj = tri_norm[i], tri_norm[j]
    if abs(dot(ni, nj)) < DOTS:
        continue
    # 共面：j 的三个顶点都在 i 的平面上
    d = dot(ni, sub(tris[j][0], tris[i][0]))
    if abs(d) > 1e-4:
        continue
    # 2D 投影求面积交叠
    u = norm(sub(tris[i][1], tris[i][0]))
    w = cross(ni, u)
    def to2(p):
        r = sub(p, tris[i][0])
        return (dot(r, u), dot(r, w))
    A, B, C = [to2(p) for p in tris[i]]
    P, Q, R = [to2(p) for p in tris[j]]
    def sgnpoly(a, b, c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    def pt_in(p, a, b, c):
        d1, d2, d3 = sgnpoly(a, b, p), sgnpoly(b, c, p), sgnpoly(c, a, p)
        E = 1e-7
        neg = (d1 < -E) or (d2 < -E) or (d3 < -E)
        pos = (d1 > E) or (d2 > E) or (d3 > E)
        return not (neg and pos)
    def area2(a, b, c): return abs(sgnpoly(a, b, c)) / 2
    a1 = area2(A, B, C); a2 = area2(P, Q, R)
    if a1 < 1e-9 or a2 < 1e-9:
        continue
    # 用 Sutherland-Hodgman 裁剪求真实交叠面积（相邻共享边的三角形交叠面积=0）
    def clip_half(poly, e0, e1):
        """保留在 e0->e1 左侧的部分"""
        if not poly: return []
        out = []
        for k in range(len(poly)):
            cur = poly[k]; nxt = poly[(k+1) % len(poly)]
            sc = sgnpoly(e0, e1, cur); sn = sgnpoly(e0, e1, nxt)
            if sc >= -1e-9: out.append(cur)
            if (sc > 1e-9 and sn < -1e-9) or (sc < -1e-9 and sn > 1e-9):
                t = sc / (sc - sn)
                out.append((cur[0]+t*(nxt[0]-cur[0]), cur[1]+t*(nxt[1]-cur[1])))
        return out
    def poly_area(poly):
        s = 0.0
        for k in range(len(poly)):
            x1, y1 = poly[k]; x2, y2 = poly[(k+1) % len(poly)]
            s += x1*y2 - x2*y1
        return abs(s)/2
    # 统一绕向为 CCW
    if sgnpoly(A, B, C) < 0: B, C = C, B
    if sgnpoly(P, Q, R) < 0: Q, R = R, Q
    tri1 = [A, B, C]; tri2 = [P, Q, R]
    res = tri2[:]
    for k in range(3):
        res = clip_half(res, tri1[k], tri1[(k+1) % 3])
        if not res: break
    inter = poly_area(res) if res else 0.0
    if inter > 1e-6:
        print("  面%d 与 面%d 共面重叠! 交叠面积 %.6f mm2 (各 %.4f / %.4f)  占比 %.1f%%"
              % (i, j, inter, a1, a2, 100*inter/min(a1, a2)))
        coplanar_dup += 1
print("  共面重叠对: %d %s" % (coplanar_dup, "(正常)" if coplanar_dup == 0 else "*** 非法 ***"))

# ---------- 2. 自相交（三角形边与其他三角形内部相交）----------
print("\n=== 2. 三角形自相交检查 ===")
def seg_tri_intersect(p0, p1, t):
    """线段 p0p1 是否与三角形 t 内部相交（不含端点共面接触）"""
    a, b, c = t
    d = p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
    e1 = sub(b, a); e2 = sub(c, a)
    h = cross(d, e2); A = dot(e1, h)
    if abs(A) < 1e-10: return False
    f = 1.0/A
    s = sub(p0, a)
    u = f * dot(s, h)
    if u < 1e-6 or u > 1-1e-6: return False
    q = cross(s, e1)
    v = f * dot(d, q)
    if v < 1e-6 or u+v > 1-1e-6: return False
    tt = f * dot(e2, q)
    return 1e-6 < tt < 1-1e-6

hits = 0
for i, ti in enumerate(tris):
    for k in range(3):
        p0, p1 = ti[k], ti[(k+1) % 3]
        for j, tj in enumerate(tris):
            if j == i: continue
            # 跳过共享顶点的三角形对
            if any(math.dist(p0, q) < 1e-6 or math.dist(p1, q) < 1e-6 for q in tj):
                continue
            if seg_tri_intersect(p0, p1, tj):
                print("  面%d 的边 与 面%d 内部相交!" % (i, j))
                hits += 1
print("  自相交数: %d %s" % (hits, "(正常)" if hits == 0 else "*** 非法 ***"))

# ---------- 3. 短边 / 退化 / 极窄三角形 ----------
print("\n=== 3. 极短边 / 极窄三角形（sliver）===")
short = []
sliver = []
for i, t in enumerate(tris):
    L = [math.dist(t[k], t[(k+1) % 3]) for k in range(3)]
    a = tri_area(t)
    if min(L) < 1e-3:
        short.append((i, min(L)))
    if a > 0:
        # 最小角
        angs = []
        for k in range(3):
            v1 = sub(t[k], t[(k+1) % 3]); v2 = sub(t[k], t[(k+2) % 3])
            c = dot(v1, v2)/(math.sqrt(dot(v1, v1))*math.sqrt(dot(v2, v2)))
            angs.append(math.degrees(math.acos(max(-1, min(1, c)))))
        if min(angs) < 1.0:
            sliver.append((i, min(angs)))
print("  边长 < 0.001mm 的三角形: %d" % len(short))
for i, L in short[:10]: print("     面%d 最短边 %.6f mm" % (i, L))
print("  最小内角 < 1° 的细长三角形: %d" % len(sliver))
for i, a in sliver[:10]: print("     面%d 最小角 %.3f°" % (i, a))

# ---------- 4. 近重合顶点（未焊接的裂缝）----------
print("\n=== 4. 近重合顶点（距离 0 但仍分离的顶点 = 裂缝隐患）===")
pts = sorted(set((round(p[0], 6), round(p[1], 6), round(p[2], 6)) for t in tris for p in t))
print("  唯一顶点数: %d" % len(pts))
close = 0
for i in range(len(pts)):
    for j in range(i+1, len(pts)):
        d = math.dist(pts[i], pts[j])
        if 1e-9 < d < 1e-3:
            print("    %.6f mm: %s <-> %s" % (d, pts[i], pts[j]))
            close += 1
print("  近重合顶点对: %d %s" % (close, "(正常)" if close == 0 else "*** 可能存在未缝合缝隙 ***"))

# ---------- 5. 壳体数（连通分量）----------
print("\n=== 5. 连通壳体数 ===")
Q = 1e-4
def vk(p): return (round(p[0]/Q), round(p[1]/Q), round(p[2]/Q))
adj = defaultdict(set)
for t in tris:
    a, b, c = vk(t[0]), vk(t[1]), vk(t[2])
    for u, v in ((a, b), (b, c), (c, a)):
        adj[u].add(v); adj[v].add(u)
seen = set(); comps = []
for s in adj:
    if s in seen: continue
    stack = [s]; comp = []
    seen.add(s)
    while stack:
        x = stack.pop(); comp.append(x)
        for y in adj[x]:
            if y not in seen:
                seen.add(y); stack.append(y)
    comps.append(comp)
print("  连通分量数: %d" % len(comps))
for i, c in enumerate(comps):
    cnt = sum(1 for t in tris if vk(t[0]) in c)
    print("     壳%d: %d 个顶点, %d 个三角形" % (i+1, len(c), cnt))
print("  %s" % ("单个封闭实体 (正常)" if len(comps) == 1 else "*** 多个独立实体! ideaMaker 可能报非法边 ***"))

# ---------- 6. 欧拉特征 ----------
V = len(adj)
E = sum(len(v) for v in adj.values()) // 2
F = len(tris)
chi = V - E + F
g = (2 - chi) // 2
print("\n=== 6. 拓扑欧拉特征 ===")
print("  V=%d  E=%d  F=%d  chi=V-E+F=%d" % (V, E, F, chi))
print("  亏格(洞数) genus = %d  -> %s" % (g, "球/带孔实体，封闭合法" if chi in (2, 0, -2) else "*** 拓扑异常 ***"))
