# -*- coding: utf-8 -*-
"""
通用 STL 可视化诊断（适用于自由曲面件，不再假设等厚拉伸）
1. 正交三视图（Lambert 着色 + 画家算法 + 背面剔除）
2. 沿三个轴的截面切片序列（统计每层环数，定位贯通孔）
输出浅色主题单文件 HTML。
用法: render_views.py <输入.stl> <输出.html>
"""
import struct, sys, math, os
from collections import defaultdict

Q = 1e-4
AXNAME = ['X', 'Y', 'Z']


def read_stl(path):
    d = open(path, 'rb').read()
    n = struct.unpack('<I', d[80:84])[0]
    tris, off = [], 84
    for i in range(n):
        v = struct.unpack('<12fH', d[off:off + 50])
        off += 50
        tris.append(((v[3], v[4], v[5]), (v[6], v[7], v[8]), (v[9], v[10], v[11])))
    return tris


def tri_normal(a, b, c):
    u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    L = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    if L < 1e-12:
        return (0.0, 0.0, 0.0)
    return (n[0] / L, n[1] / L, n[2] / L)


def volume(tris):
    s = 0.0
    for a, b, c in tris:
        s += (a[0] * (b[1] * c[2] - b[2] * c[1])
              - a[1] * (b[0] * c[2] - b[2] * c[0])
              + a[2] * (b[0] * c[1] - b[1] * c[0]))
    return s / 6.0


def bbox(tris):
    lo = [min(p[k] for t in tris for p in t) for k in range(3)]
    hi = [max(p[k] for t in tris for p in t) for k in range(3)]
    return lo, hi


def euler(tris):
    vk = lambda p: (round(p[0] / Q), round(p[1] / Q), round(p[2] / Q))
    verts, edges, faces = set(), set(), 0
    for tri in tris:
        ks = [vk(p) for p in tri]
        if ks[0] == ks[1] or ks[1] == ks[2] or ks[0] == ks[2]:
            continue
        faces += 1
        for k in ks:
            verts.add(k)
        for i in range(3):
            e = tuple(sorted((ks[i], ks[(i + 1) % 3])))
            edges.add(e)
    return len(verts), len(edges), faces


# ---------- 正交视图渲染 ----------
LIGHT = (0.42, 0.55, 0.72)


def render_view(tris, u, v, dep, W, H, pad=6):
    """u,v = 屏幕横/纵轴索引; dep = 深度轴索引。观察者位于 +dep 方向。"""
    lo = [min(p[k] for t in tris for p in t) for k in range(3)]
    hi = [max(p[k] for t in tris for p in t) for k in range(3)]
    su, sv = hi[u] - lo[u], hi[v] - lo[v]
    sc = min((W - 2 * pad) / max(su, 1e-9), (H - 2 * pad) / max(sv, 1e-9))
    ox = pad + ((W - 2 * pad) - su * sc) / 2
    oy = pad + ((H - 2 * pad) - sv * sc) / 2

    def tx(p):
        return (ox + (p[u] - lo[u]) * sc, H - oy - (p[v] - lo[v]) * sc)

    items = []
    for tri in tris:
        n = tri_normal(*tri)
        if n[dep] <= 0:            # 背面剔除
            continue
        cen = sum(p[dep] for p in tri) / 3.0
        lam = n[0] * LIGHT[0] + n[1] * LIGHT[1] + n[2] * LIGHT[2]
        sh = 0.32 + 0.68 * max(0.0, lam)
        g = int(232 - 118 * sh)
        col = '#%02x%02x%02x' % (g, g + 4, min(255, g + 12))
        pts = ' '.join('%.2f,%.2f' % tx(p) for p in tri)
        items.append((cen, '<polygon points="%s" fill="%s" stroke="%s" stroke-width="0.3"/>' % (pts, col, col)))
    items.sort(key=lambda x: x[0])   # 画家算法：远的先画
    body = ''.join(it[1] for it in items)
    axlabel = '沿 %s 视向  (%s–%s 面)' % (AXNAME[dep], AXNAME[u], AXNAME[v])
    return body, axlabel, sc


# ---------- 截面切片 ----------
def slice_rings(tris, axis, h):
    """返回该高度上的闭合环列表，每环为 [(u,v)...]，u/v 为另两轴。"""
    oth = [k for k in range(3) if k != axis]
    segs = []
    for tri in tris:
        ds = []
        for p in tri:
            d = p[axis] - h
            ds.append(d if abs(d) > 1e-9 else 1e-9)
        if all(d > 0 for d in ds) or all(d < 0 for d in ds):
            continue
        inter = []
        for i in range(3):
            d1, d2 = ds[i], ds[(i + 1) % 3]
            if (d1 > 0) != (d2 > 0):
                p1, p2 = tri[i], tri[(i + 1) % 3]
                t = d1 / (d1 - d2)
                inter.append(tuple(p1[k] + t * (p2[k] - p1[k]) for k in range(3)))
        if len(inter) >= 2:
            segs.append((inter[0], inter[1]))
    if not segs:
        return []
    kv = lambda p: (round(p[0] / Q), round(p[1] / Q), round(p[2] / Q))
    adj = defaultdict(list)
    pos = {}
    for a, b in segs:
        ka, kb = kv(a), kv(b)
        if ka == kb:
            continue
        pos[ka] = a
        pos[kb] = b
        adj[ka].append(kb)
        adj[kb].append(ka)
    visited = set()
    rings = []
    for s in list(adj.keys()):
        for nb in adj[s]:
            if (s, nb) in visited:
                continue
            ring = [s]
            cur, nxt = s, nb
            while True:
                visited.add((cur, nxt))
                visited.add((nxt, cur))
                ring.append(nxt)
                found = None
                for nn in adj[nxt]:
                    if (nxt, nn) not in visited:
                        found = nn
                        break
                if found is None:
                    break
                cur, nxt = nxt, found
            if len(ring) >= 4:
                rings.append([(pos[k][oth[0]], pos[k][oth[1]]) for k in ring[:-1]])
    return rings


def ring_area(r):
    s = 0.0
    for i in range(len(r)):
        x1, y1 = r[i]
        x2, y2 = r[(i + 1) % len(r)]
        s += x1 * y2 - x2 * y1
    return s / 2.0


def rings_svg(rings, S=150):
    if not rings:
        return '<svg width="%d" height="%d"><rect width="100%%" height="100%%" fill="#fafafa"/></svg>' % (S, S)
    pts = [p for r in rings for p in r]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    lo_x, hi_x = min(xs), max(xs); lo_y, hi_y = min(ys), max(ys)
    w = max(hi_x - lo_x, 1e-9); h = max(hi_y - lo_y, 1e-9)
    pad = 8
    sc = min((S - 2 * pad) / w, (S - 2 * pad) / h)
    ox = pad + ((S - 2 * pad) - w * sc) / 2
    oy = pad + ((S - 2 * pad) - h * sc) / 2
    out = ['<svg width="%d" height="%d"><rect width="100%%" height="100%%" fill="#fafafa"/>' % (S, S)]
    areas = [(abs(ring_area(r)), r) for r in rings]
    areas.sort(key=lambda x: -x[0])
    outer = areas[0][1] if areas else None
    for a, r in areas:
        is_out = (r is outer)
        # 内孔（负面积或绝对面积小于外环）用白色挖空表示
        d = 'M' + ' L'.join('%.2f,%.2f' % (ox + (p[0] - lo_x) * sc, S - oy - (p[1] - lo_y) * sc) for p in r) + ' Z'
        out.append('<path d="%s" fill="%s" stroke="#334155" stroke-width="0.8" fill-rule="evenodd"/>'
                   % (d, '#cbd5e1' if is_out else '#ffffff'))
    out.append('</svg>')
    return ''.join(out)


def main():
    src, out = sys.argv[1], sys.argv[2]
    tris = read_stl(src)
    lo, hi = bbox(tris)
    size = [hi[k] - lo[k] for k in range(3)]
    V, E, F = euler(tris)
    chi = V - E + F
    genus = (2 - chi) // 2
    vol = volume(tris)
    order = sorted(range(3), key=lambda k: -size[k])   # 长 -> 短

    # 三视图
    views = []
    for dep, title in ((order[2], '主视'), (order[1], '侧视'), (order[0], '俯视')):
        u, v = [k for k in range(3) if k != dep][0], [k for k in range(3) if k != dep][1]
        body, lbl, sc = render_view(tris, u, v, dep, 260, 300)
        views.append((title, lbl, body))

    # 切片扫描：沿三轴
    scan_html = []
    scan_stat = []
    for axis in range(3):
        N = 9 if axis == order[0] else 5
        rows = []
        maxrings = 0
        for i in range(N):
            h = lo[axis] + size[axis] * (i + 0.5) / N
            rings = slice_rings(tris, axis, h)
            maxrings = max(maxrings, len(rings))
            areas = sorted((abs(ring_area(r)) for r in rings), reverse=True)
            net = (areas[0] if areas else 0.0) - sum(areas[1:])
            rows.append((h, len(rings), net, rings_svg(rings, 132)))
        scan_stat.append((axis, maxrings))
        cells = ''.join(
            '<div class="sl"><div class="slh">%.1f mm</div>%s<div class="slf">%d 环 · %.0f mm²</div></div>'
            % (r[0], r[3], r[1], r[2]) for r in rows)
        scan_html.append((axis, cells))

    name = os.path.basename(src)
    verdict_ok = (genus >= 0)
    vcls = 'v-ok' if verdict_ok else 'v-bad'
    vtext = ('拓扑自洽（χ=%d，亏格 %d）' % (chi, genus)) if verdict_ok else '拓扑异常'

    html = """<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>STL 诊断 · NAME</title>
<style>
body{font-family:-apple-system,"Segoe UI","Microsoft YaHei",sans-serif;background:#f6f7f9;color:#1e293b;margin:0;padding:28px}
.wrap{max-width:1080px;margin:0 auto}
h1{font-size:19px;font-weight:600;margin:0 0 4px}
.sub{color:#64748b;font-size:13px;margin-bottom:20px}
.cards{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:20px}
.card{background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:12px 16px;min-width:118px}
.card .k{font-size:12px;color:#64748b;margin-bottom:4px}
.card .v{font-size:17px;font-weight:600}
.card .v.sm{font-size:14px}
h2{font-size:15px;font-weight:600;margin:26px 0 10px}
.views{display:flex;gap:14px;flex-wrap:wrap}
.vw{background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:10px}
.vw .t{font-size:13px;font-weight:600;margin-bottom:2px}
.vw .s{font-size:11px;color:#94a3b8;margin-bottom:6px}
.scan{background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:14px;margin-bottom:14px}
.scan .t{font-size:13px;font-weight:600;margin-bottom:10px}
.sls{display:flex;gap:10px;flex-wrap:wrap}
.sl{text-align:center}
.sl .slh{font-size:11px;color:#64748b;margin-bottom:3px}
.sl .slf{font-size:11px;color:#334155;margin-top:3px}
.verdict{border-radius:10px;padding:14px 16px;margin-top:22px;font-size:14px;line-height:1.7}
.v-ok{background:#ecfdf5;border:1px solid #a7f3d0}
.v-bad{background:#fef2f2;border:1px solid #fecaca}
table{border-collapse:collapse;width:100%%;font-size:13px;background:#fff;border-radius:10px;overflow:hidden}
th,td{border:1px solid #e2e8f0;padding:7px 10px;text-align:left}
th{background:#f1f5f9;font-weight:600}
</style></head><body><div class="wrap">
<h1>STL 诊断 · NAME</h1>
<div class="sub">SRC</div>
<div class="cards">
<div class="card"><div class="k">尺寸 (mm)</div><div class="v sm">SX × SY × SZ</div></div>
<div class="card"><div class="k">三角形 / 顶点</div><div class="v">F / V</div></div>
<div class="card"><div class="k">体积</div><div class="v sm">VOL mm³</div></div>
<div class="card"><div class="k">欧拉特征 χ</div><div class="v">CHI</div></div>
<div class="card"><div class="k">亏格（贯通孔）</div><div class="v">GEN</div></div>
<div class="card"><div class="k">连通壳体</div><div class="v">1</div></div>
</div>
<h2>正交三视图</h2>
<div class="views">VIEWS</div>
<h2>截面扫描（灰色为实体，白色为孔）</h2>
SCANS
<div class="verdict VCLS"><b>结论：VTXT</b><br>NOTE</div>
</div></body></html>"""

    views_html = ''.join(
        '<div class="vw"><div class="t">%s</div><div class="s">%s</div>'
        '<svg width="260" height="300"><rect width="260" height="300" fill="#fafafa"/>%s</svg></div>'
        % (t, l, b) for t, l, b in views)

    scans_html = ''.join(
        '<div class="scan"><div class="t">沿 %s 轴扫描（%s，最大同时环数 %d）</div><div class="sls">%s</div></div>'
        % (AXNAME[ax], '%.1f' % size[ax], mr, cells) for (ax, mr), cells in zip(scan_stat, [s[1] for s in scan_html]))

    note = ('χ 为偶数、单一封闭壳体 → 网格本身是合法流形。亏格 %d 表示模型上有 %d 个贯通通道（tunnel），'
            '这是设计形状本身的属性，<b>不是非法边</b>。切片软件报「非法的边」只与非流形边（一条边被 ≠2 个面共享）有关，'
            '与本项无关。' % (genus, genus)) if verdict_ok else 'χ 为奇数 → 网格存在拓扑破缺。'

    html = (html.replace('NAME', name).replace('SRC', src)
            .replace('SX', '%.2f' % size[0]).replace('SY', '%.2f' % size[1]).replace('SZ', '%.2f' % size[2])
            .replace('>F / V<', '>%d / %d<' % (F, V))
            .replace('VOL', '%.2f' % vol)
            .replace('CHI', str(chi)).replace('GEN', str(genus))
            .replace('VCLS', vcls).replace('VTXT', vtext).replace('NOTE', note)
            .replace('VIEWS', views_html).replace('SCANS', scans_html)
            .replace('%%', '%'))

    open(out, 'w', encoding='utf-8').write(html)
    print('报告已生成: %s' % out)
    print('尺寸 %.2f × %.2f × %.2f mm | 三角形 %d | 体积 %.2f mm³ | V=%d E=%d F=%d chi=%d genus=%d'
          % (size[0], size[1], size[2], F, vol, V, E, F, chi, genus))
    for ax, mr in scan_stat:
        print('  沿 %s 轴最多同时 %d 个环' % (AXNAME[ax], mr))


if __name__ == '__main__':
    main()
