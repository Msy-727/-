# -*- coding: utf-8 -*-
"""
STL 终端 ASCII 正交预览（纯标准库，无需 PIL/matplotlib）
用法: ascii_view.py <输入.stl> [字符宽度]
输出沿 X / Y / Z 三个视向的正交投影，z-buffer 光栅化 + Lambert 明暗。
"""
import struct, sys, math

RAMP = ' .:-=+*#%@'
LIGHT = (0.42, 0.55, 0.72)


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


def render(tris, dep, W=78):
    u, v = [k for k in range(3) if k != dep]
    lo = [min(p[k] for t in tris for p in t) for k in range(3)]
    hi = [max(p[k] for t in tris for p in t) for k in range(3)]
    su, sv = hi[u] - lo[u], hi[v] - lo[v]
    if su < 1e-9 or sv < 1e-9:
        return ['(该方向退化)']
    H = max(6, int(round(W * (sv / su) * 0.5)))
    zbuf = [[-1e30] * W for _ in range(H)]
    cbuf = [[' '] * W for _ in range(H)]
    sx = (W - 1) / su
    sy = (H - 1) / sv
    for tri in tris:
        n = tri_normal(*tri)
        if n[dep] <= 0:
            continue
        lam = n[0] * LIGHT[0] + n[1] * LIGHT[1] + n[2] * LIGHT[2]
        sh = 0.25 + 0.75 * max(0.0, lam)
        ch = RAMP[min(len(RAMP) - 1, int(sh * (len(RAMP) - 1) + 0.5))]
        A = ((tri[0][u] - lo[u]) * sx, (tri[0][v] - lo[v]) * sy, tri[0][dep])
        B = ((tri[1][u] - lo[u]) * sx, (tri[1][v] - lo[v]) * sy, tri[1][dep])
        C = ((tri[2][u] - lo[u]) * sx, (tri[2][v] - lo[v]) * sy, tri[2][dep])
        minx = max(0, int(math.floor(min(A[0], B[0], C[0]))))
        maxx = min(W - 1, int(math.ceil(max(A[0], B[0], C[0]))))
        miny = max(0, int(math.floor(min(A[1], B[1], C[1]))))
        maxy = min(H - 1, int(math.ceil(max(A[1], B[1], C[1]))))
        d = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1])
        if abs(d) < 1e-12:
            continue
        for py in range(miny, maxy + 1):
            for px in range(minx, maxx + 1):
                x, y = px + 0.5, py + 0.5
                w0 = ((B[0] - A[0]) * (y - A[1]) - (x - A[0]) * (B[1] - A[1])) / d
                w1 = ((x - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (y - A[1])) / d
                if w0 < -1e-9 or w1 < -1e-9:
                    continue
                w2 = 1.0 - w0 - w1
                if w2 < -1e-9:
                    continue
                z = w2 * A[2] + w1 * B[2] + w0 * C[2]
                ry = H - 1 - py
                if z > zbuf[ry][px]:
                    zbuf[ry][px] = z
                    cbuf[ry][px] = ch
    return [''.join(r) for r in cbuf]


def main():
    src = sys.argv[1]
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 78
    tris = read_stl(src)
    lo = [min(p[k] for t in tris for p in t) for k in range(3)]
    hi = [max(p[k] for t in tris for p in t) for k in range(3)]
    size = [hi[k] - lo[k] for k in range(3)]
    print('尺寸 X %.2f  Y %.2f  Z %.2f  | 三角形 %d' % (size[0], size[1], size[2], len(tris)))
    order = sorted(range(3), key=lambda k: -size[k])
    for dep in order:
        an = 'XYZ'[dep]
        o1, o2 = [k for k in range(3) if k != dep]
        print('\n=== 沿 %s 视向  (横轴 %s = %.1f mm, 纵轴 %s = %.1f mm) ==='
              % (an, 'XYZ'[o1], size[o1], 'XYZ'[o2], size[o2]))
        for line in render(tris, dep, W):
            print(line)


if __name__ == '__main__':
    main()
