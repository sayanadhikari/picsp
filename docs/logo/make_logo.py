"""Generates the PICSP3D logo (SVG + PNG): an isometric mesh cube threaded by a
gyrating particle orbit, with the wordmark. Run: python3 make_logo.py"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY, INK, CYAN, ORANGE, GRID = "#0f1b2d", "#0f1b2d", "#19b5d8", "#f2662f", "#7d8ea6"


def proj(p, s, ox, oy):
    x, y, z = p
    c, d = math.cos(math.radians(30)), math.sin(math.radians(30))
    return ox + s * (x - y) * c, oy + s * ((x + y) * d - z)


def emblem(s, ox, oy, dark):
    out = []
    edge = "#e8eef6" if dark else NAVY
    P = lambda p: proj(p, s, ox, oy)
    line = lambda a, b, col, w, op=1.0, dash="": out.append(
        f'<line x1="{P(a)[0]:.2f}" y1="{P(a)[1]:.2f}" x2="{P(b)[0]:.2f}" y2="{P(b)[1]:.2f}" stroke="{col}" '
        f'stroke-width="{w}" stroke-opacity="{op}" stroke-linecap="round"{dash}/>')
    # three visible faces with a light fill
    faces = [((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)), ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),
             ((0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1))]
    shades = (0.10, 0.05, 0.07) if dark else (0.07, 0.035, 0.05)
    for f, sh in zip(faces, shades):
        pts = " ".join(f"{P(q)[0]:.2f},{P(q)[1]:.2f}" for q in f)
        out.append(f'<polygon points="{pts}" fill="{CYAN}" fill-opacity="{sh}"/>')
    # hidden edges
    for a, b in [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1))]:
        line(a, b, GRID, s * 0.006, 0.55, ' stroke-dasharray="3 4"')
    # mesh lines on the visible faces
    n = 4
    for k in range(1, n):
        t = k / n
        line((t, 0, 1), (t, 1, 1), GRID, s * 0.004, 0.6)
        line((0, t, 1), (1, t, 1), GRID, s * 0.004, 0.6)
        line((1, t, 0), (1, t, 1), GRID, s * 0.004, 0.6)
        line((1, 0, t), (1, 1, t), GRID, s * 0.004, 0.6)
        line((t, 1, 0), (t, 1, 1), GRID, s * 0.004, 0.6)
        line((0, 1, t), (1, 1, t), GRID, s * 0.004, 0.6)
    # gyro orbit: a helix along z, slightly tilted, from below the cube to above it
    pts = []
    N = 420
    for i in range(N + 1):
        u = i / N
        z = -0.18 + 1.36 * u
        ph = 2 * math.pi * 4.25 * u
        r = 0.30 - 0.06 * math.sin(math.pi * u)  # mirror-like neck in the middle
        x, y = 0.5 + r * math.cos(ph) + 0.05 * (u - 0.5), 0.5 + r * math.sin(ph) - 0.05 * (u - 0.5)
        depth = math.cos(ph - math.radians(45))  # > 0 towards the viewer
        pts.append(((x, y, z), depth, u))
    back, front = [], []
    for i in range(N):
        a, b = pts[i], pts[i + 1]
        (back if a[1] < 0 else front).append((a, b))

    def seg(a, b, op, w):
        u = a[2]
        col = mix(CYAN, ORANGE, u)
        line(a[0], b[0], col, w, op)

    for a, b in back:
        seg(a, b, 0.35, s * 0.016)
    # front cube edges drawn between the back and the front part of the orbit
    for a, b in [((1, 0, 0), (1, 1, 0)), ((1, 1, 0), (0, 1, 0)), ((1, 1, 0), (1, 1, 1)), ((1, 0, 0), (1, 0, 1)),
                 ((0, 1, 0), (0, 1, 1)), ((0, 0, 1), (1, 0, 1)), ((1, 0, 1), (1, 1, 1)), ((1, 1, 1), (0, 1, 1)),
                 ((0, 1, 1), (0, 0, 1))]:
        line(a, b, edge, s * 0.014, 1.0)
    for a, b in front:
        seg(a, b, 1.0, s * 0.024)
    # guiding-centre / field line
    line((0.5, 0.5, -0.3), (0.5, 0.5, 1.32), edge, s * 0.006, 0.45, ' stroke-dasharray="2 5"')
    # the particle at the top of its orbit, with a glow
    hx, hy = P(pts[-1][0])
    out.append(f'<circle cx="{hx:.2f}" cy="{hy:.2f}" r="{s * 0.075:.2f}" fill="{ORANGE}" fill-opacity="0.22"/>')
    out.append(f'<circle cx="{hx:.2f}" cy="{hy:.2f}" r="{s * 0.042:.2f}" fill="{ORANGE}"/>')
    # a few macro-particles on the mesh
    for q, c in [((0.25, 0.75, 1), CYAN), ((0.75, 0.25, 1), ORANGE), ((1, 0.75, 0.25), CYAN), ((0.5, 1, 0.75), CYAN),
                 ((1, 0.25, 0.75), ORANGE)]:
        x, y = P(q)
        out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{s * 0.022:.2f}" fill="{c}"/>')
    return "\n".join(out)


def mix(a, b, t):
    ca = [int(a[i:i + 2], 16) for i in (1, 3, 5)]
    cb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(ca, cb))


def svg(w, h, body, bg=None):
    rect = f'<rect width="{w}" height="{h}" rx="{min(w, h) * 0.08:.0f}" fill="{bg}"/>' if bg else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
            f'{rect}{body}</svg>')


FONT = "font-family=\"'Montserrat','Inter','Helvetica Neue',Arial,sans-serif\""


def wordmark(x, y, size, dark):
    ink = "#f4f7fb" if dark else INK
    sub = "#9fb0c6" if dark else "#51627a"
    return (f'<text x="{x}" y="{y}" {FONT} font-size="{size}" font-weight="800" letter-spacing="{size * 0.02:.1f}" '
            f'fill="{ink}">PICSP<tspan fill="{ORANGE}">3D</tspan></text>'
            f'<text x="{x + size * 0.04:.1f}" y="{y + size * 0.42:.1f}" {FONT} font-size="{size * 0.128:.1f}" '
            f'font-weight="500" letter-spacing="{size * 0.03:.1f}" fill="{sub}">3D3V PARTICLE-IN-CELL FOR PLASMA-OBJECT INTERACTION</text>')


def main():
    files = {}
    files["picsp3d_logo.svg"] = svg(1200, 360, emblem(125, 150, 180, False) + wordmark(300, 200, 150, False))
    files["picsp3d_logo_dark.svg"] = svg(1200, 360, emblem(125, 150, 180, True) + wordmark(300, 200, 150, True), NAVY)
    # transparent version with light lettering for dark pages (GitHub dark mode)
    files["picsp3d_logo_on_dark.svg"] = svg(1200, 360, emblem(125, 150, 180, True) + wordmark(300, 200, 150, True))
    files["picsp3d_icon.svg"] = svg(512, 512, emblem(170, 256, 262, True), NAVY)
    import cairosvg
    for name, text in files.items():
        p = os.path.join(HERE, name)
        open(p, "w").write(text)
        cairosvg.svg2png(bytestring=text.encode(), write_to=p.replace(".svg", ".png"), output_width=
                         1200 if "icon" not in name else 512)
        print("wrote", name)


if __name__ == "__main__":
    main()
