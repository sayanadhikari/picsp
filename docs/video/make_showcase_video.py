"""PICSP3D showcase video (square 1080 x 1080, H.264) for social media.

Built from the feature-example figures. For sharp animations render the example frames at high
resolution first (lossless PNG frames, documentation media untouched):

    cd examples/features
    export PICSP_MEDIA_DIR=/tmp/hires PICSP_FRAMES=1 PICSP_ANIM_DPI=200
    for e in 02_two_stream_instability 04_magnetic_mirror 08_sphere_wake 12_spacecraft_STL; do
        python3 run_example.py $e --analyze-only
    done
    cd ../../docs/video
    python3 make_showcase_video.py --frames /tmp/hires -o picsp3d_showcase.mp4

Without --frames the GIFs in docs/media are used.
"""
import argparse
import glob
import math
import os
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MEDIA = os.path.join(ROOT, "docs", "media")
LOGO = os.path.join(ROOT, "docs", "logo", "picsp3d_logo_on_dark.png")

W = H = 1080
FPS = 30
BG0, BG1 = (9, 16, 30), (18, 34, 58)
INK, SUB, MUTED = (244, 247, 251), (170, 186, 206), (110, 128, 152)
ORANGE, CYAN = (242, 102, 47), (25, 181, 216)


def font(weight, size):
    for p in (f"/usr/share/fonts/opentype/inter/Inter-{weight}.otf", f"/usr/share/fonts/truetype/inter/Inter-{weight}.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if weight in ("Bold", "ExtraBold")
              else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


F_KICK = font("SemiBold", 22)
F_HEAD = font("Bold", 60)
F_SUB = font("Regular", 30)
F_CHIP = font("SemiBold", 34)
F_MONO = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 24) \
    if os.path.exists("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf") else F_SUB


def ease(t):
    t = min(max(t, 0.0), 1.0)
    return t * t * (3 - 2 * t)


# ---------------------------------------------------------------------------------------
# background: dark gradient with a slowly drifting PIC mesh and a few particles
# ---------------------------------------------------------------------------------------
_grad = None


def background(t):
    global _grad
    if _grad is None:
        y = np.linspace(0, 1, H)[:, None, None]
        x = np.linspace(0, 1, W)[None, :, None]
        a = 0.65 * y + 0.35 * x
        _grad = (np.array(BG0) * (1 - a) + np.array(BG1) * a).astype(np.uint8)
    im = Image.fromarray(np.broadcast_to(_grad, (H, W, 3)).copy())
    d = ImageDraw.Draw(im, "RGBA")
    step = 60
    off = (t * 6) % step
    for k in range(-1, W // step + 2):
        x = k * step + off
        d.line([(x, 0), (x, H)], fill=(120, 150, 190, 14), width=1)
    for k in range(-1, H // step + 2):
        y = k * step + off * 0.5
        d.line([(0, y), (W, y)], fill=(120, 150, 190, 14), width=1)
    rng = np.random.default_rng(3)
    for _ in range(40):
        x0, y0, vx, vy, r = rng.uniform(0, W), rng.uniform(0, H), rng.normal(0, 14), rng.normal(0, 14), rng.uniform(1.2, 2.6)
        c = ORANGE if rng.random() < 0.35 else CYAN
        x, y = (x0 + vx * t) % W, (y0 + vy * t) % H
        d.ellipse([x - r, y - r, x + r, y + r], fill=c + (70,))
    return im


# ---------------------------------------------------------------------------------------
# media helpers
# ---------------------------------------------------------------------------------------
def load_frames(name, frames_dir):
    """Frames of an animation: PNG frames from --frames if present, else the GIF in docs/media."""
    if frames_dir:
        files = sorted(glob.glob(os.path.join(frames_dir, f"{name}_frames", "*.png")))
        if files:
            return [Image.open(f).convert("RGB") for f in files]
    gif = Image.open(os.path.join(MEDIA, f"{name}.gif"))
    out = []
    for i in range(getattr(gif, "n_frames", 1)):
        gif.seek(i)
        out.append(gif.convert("RGB"))
    return out


def still(name):
    return Image.open(os.path.join(MEDIA, name)).convert("RGB")


def fit(img, box_w, box_h):
    s = min(box_w / img.width, box_h / img.height)
    return img.resize((max(1, int(img.width * s)), max(1, int(img.height * s))), Image.LANCZOS)


def card(canvas, img, cx, cy, box_w, box_h, alpha=1.0, light=True, zoom=1.0, radius=22):
    """Paste `img` scaled into a rounded card centred at (cx, cy)."""
    im = fit(img, box_w, box_h)
    if zoom != 1.0:  # slow push-in: crop the centre
        zw, zh = int(im.width / zoom), int(im.height / zoom)
        x0, y0 = (im.width - zw) // 2, (im.height - zh) // 2
        im = im.crop((x0, y0, x0 + zw, y0 + zh)).resize(im.size, Image.LANCZOS)
    pad = 18 if light else 0
    cw, ch = im.width + 2 * pad, im.height + 2 * pad
    base = Image.new("RGBA", (cw, ch), (0, 0, 0, 0))
    mask = Image.new("L", (cw, ch), 0)
    ImageDraw.Draw(mask).rounded_rectangle([0, 0, cw - 1, ch - 1], radius, fill=255)
    bgc = (250, 250, 248, 255) if light else (0, 0, 0, 255)
    base.paste(Image.new("RGBA", (cw, ch), bgc), (0, 0), mask)
    base.paste(im, (pad, pad))
    base.putalpha(Image.fromarray((np.asarray(mask, float) * alpha).astype(np.uint8)))
    # soft shadow
    sh = Image.new("RGBA", (cw + 60, ch + 60), (0, 0, 0, 0))
    ImageDraw.Draw(sh).rounded_rectangle([30, 38, cw + 30, ch + 38], radius, fill=(0, 0, 0, int(120 * alpha)))
    sh = sh.filter(ImageFilter.GaussianBlur(18))
    x, y = int(cx - cw / 2), int(cy - ch / 2)
    canvas.alpha_composite(sh, (x - 30, y - 30))
    canvas.alpha_composite(base, (x, y))
    if not light:
        ImageDraw.Draw(canvas).rounded_rectangle([x, y, x + cw - 1, y + ch - 1], radius,
                                                 outline=(120, 150, 190, int(90 * alpha)), width=2)


def text(d, xy, s, f, fill, alpha=1.0, anchor="la"):
    d.text(xy, s, font=f, fill=fill + (int(255 * alpha),), anchor=anchor)


def wrap(s, f, width, d):
    words, lines, cur = s.split(), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if d.textlength(t, font=f) <= width:
            cur = t
        else:
            lines.append(cur)
            cur = w
    return lines + ([cur] if cur else [])


def titles(canvas, local, dur, kicker, head, sub):
    """Kicker, headline and subline with a staggered slide-up fade."""
    d = ImageDraw.Draw(canvas, "RGBA")
    out = 1 - ease((local - (dur - 0.45)) / 0.45)
    a1, a2, a3 = (ease((local - k) / 0.5) * out for k in (0.0, 0.12, 0.28))
    y = 92
    d.rounded_rectangle([72, y + 2, 78, y + 26], 3, fill=ORANGE + (int(255 * a1),))
    text(d, (92, y + int(12 * (1 - a1))), kicker.upper(), F_KICK, ORANGE, a1)
    y += 46
    for line in wrap(head, F_HEAD, W - 144, d):
        text(d, (72, y + int(16 * (1 - a2))), line, F_HEAD, INK, a2)
        y += 70
    y += 6
    for line in wrap(sub, F_SUB, W - 144, d):
        text(d, (72, y + int(14 * (1 - a3))), line, F_SUB, SUB, a3)
        y += 40
    return y


def progress(canvas, t, total):
    d = ImageDraw.Draw(canvas, "RGBA")
    d.rectangle([0, H - 6, W * t / total, H], fill=ORANGE + (200,))
    text(d, (W - 72, H - 44), "PICSP3D", F_KICK, MUTED, 0.9, anchor="ra")


# ---------------------------------------------------------------------------------------
# scenes: each is (duration, draw(canvas, local_time, dur))
# ---------------------------------------------------------------------------------------
def build_scenes(frames_dir):
    logo = Image.open(LOGO).convert("RGBA")
    S = []

    def intro(c, t, dur):
        a = ease(t / 0.9) * (1 - ease((t - (dur - 0.5)) / 0.5))
        s = 0.92 + 0.08 * ease(t / 1.4)
        lg = logo.resize((int(940 * s), int(940 * s * logo.height / logo.width)), Image.LANCZOS)
        lg.putalpha(Image.fromarray((np.asarray(lg.getchannel("A"), float) * a).astype(np.uint8)))
        c.alpha_composite(lg, ((W - lg.width) // 2, 330))
        d = ImageDraw.Draw(c, "RGBA")
        b = ease((t - 0.8) / 0.6) * (1 - ease((t - (dur - 0.5)) / 0.5))
        text(d, (W // 2, 690), "Open-source plasma simulation", F_HEAD, INK, b, anchor="ma")
        text(d, (W // 2, 775), "Spacecraft charging · sheaths and wakes · probes · magnetised plasmas",
             font("Regular", 26), SUB, b, anchor="ma")
    S.append((4.2, intro))

    def still_scene(kicker, head, sub, img, light, box_h=560, zoom_to=1.06):
        def f(c, t, dur):
            y = titles(c, t, dur, kicker, head, sub)
            a = ease((t - 0.25) / 0.6) * (1 - ease((t - (dur - 0.45)) / 0.45))
            cy = y + 40 + box_h / 2
            card(c, img, W / 2, min(cy, H - 70 - box_h / 2), W - 120, box_h, a, light, 1 + (zoom_to - 1) * t / dur)
        return f

    def anim_scene(kicker, head, sub, frames, fps, light=True, box_h=560):
        def f(c, t, dur):
            y = titles(c, t, dur, kicker, head, sub)
            a = ease((t - 0.25) / 0.6) * (1 - ease((t - (dur - 0.45)) / 0.45))
            i = min(len(frames) - 1, int(max(0.0, t - 0.2) * fps))  # play once, hold the last frame
            cy = y + 40 + box_h / 2
            card(c, frames[i], W / 2, min(cy, H - 70 - box_h / 2), W - 120, box_h, a, light)
        return f

    S.append((4.6, still_scene("Spacecraft charging", "A CubeSat in low Earth orbit",
                               "Imported from CAD with solar panels and an isolated probe. Sunlit faces photoemit, "
                               "shaded faces collect electrons.", still("12_cubesat_render.png"), False, 600)))
    S.append((4.6, anim_scene("Plasma flow", "Ion wakes behind spacecraft",
                              "Mesothermal O⁺ flow: the wake empties and refills as ions focus downstream.",
                              load_frames("12_cubesat_wake", frames_dir), 4.0, True, 600)))
    S.append((4.4, anim_scene("Plasma flow", "Wake of a floating sphere",
                              "Ion void, focusing and a positive wake potential.",
                              load_frames("08_wake", frames_dir), 8.0, True, 520)))
    S.append((5.4, anim_scene("Kinetic physics", "Two-stream instability",
                              "Counter-streaming beams trap each other into phase-space vortices. Growth rate "
                              "0.48 vs 0.50 ω_b from theory.", load_frames("02_two_stream_phase", frames_dir), 10.0, True,
                              560)))
    S.append((4.6, anim_scene("Magnetised plasma", "A magnetic bottle from two coils",
                              "Exact current-loop fields. Trapped electrons bounce, the loss cone empties: "
                              "trapped fraction 0.841 vs 0.839.", load_frames("04_mirror_bounce", frames_dir), 10.0,
                              True, 500)))
    S.append((4.2, still_scene("Magnetised plasma", "Protons in a dipole field",
                               "Gyration, bounce between mirror points and azimuthal drift.",
                               still("05_dipole_trapping.png"), False, 600, 1.12)))
    S.append((4.2, still_scene("Surface physics", "Photoemission with self-shadowing",
                               "Ray-traced sunlight: the mast shades the panel. Secondary emission too.",
                               still("11_photoemission_shadow.png"), False, 600, 1.08)))

    def collisions(c, t, dur):
        y = titles(c, t, dur, "Collisions", "Neutral and Coulomb collisions",
                   "Monte Carlo collisions, ionisation avalanches and binary Coulomb relaxation, checked "
                   "against theory.")
        out = 1 - ease((t - (dur - 0.45)) / 0.45)
        a1, a2 = ease((t - 0.25) / 0.6) * out, ease((t - 0.9) / 0.6) * out
        card(c, still("14_ionization_avalanche.png"), W / 2, y + 190, W - 120, 300, a1)
        card(c, still("15_coulomb_relaxation.png"), W / 2, y + 520, W - 120, 300, a2)
    S.append((5.0, collisions))

    def cad(c, t, dur):
        y = titles(c, t, dur, "Geometry", "From CAD file to simulation",
                   "Units, placement and grid in one command.")
        d = ImageDraw.Draw(c, "RGBA")
        out = 1 - ease((t - (dur - 0.45)) / 0.45)
        cmd = "$ picsp-post mesh import part.stl --units mm"
        n = int(len(cmd) * ease((t - 0.3) / 1.6))
        a = ease((t - 0.2) / 0.4) * out
        d.rounded_rectangle([60, y + 10, W - 60, y + 74], 14, fill=(6, 11, 22, int(235 * a)),
                            outline=(120, 150, 190, int(110 * a)), width=2)
        text(d, (84, y + 28), cmd[:n] + ("█" if int(t * 2) % 2 == 0 and n < len(cmd) else ""), F_MONO,
             (180, 230, 200), a)
        a2 = ease((t - 1.4) / 0.6) * out
        card(c, still("16_mesh_import.png"), W / 2, y + 330, W - 120, 420, a2)
    S.append((5.2, cad))

    chips = ["MPI + OpenMP", "Multigrid Poisson", "Cut-cell boundaries", "Floating conductors",
             "STL + CSG geometry", "Static and varying B", "Photo + secondary emission", "MCC + Coulomb collisions",
             "ParaView output", "Python post-processing"]

    def features(c, t, dur):
        y = titles(c, t, dur, "Under the hood", "Built for HPC", "Fast 3D3V electrostatic PIC in modern C++.")
        d = ImageDraw.Draw(c, "RGBA")
        out = 1 - ease((t - (dur - 0.45)) / 0.45)
        x, yy = 72, y + 60
        for k, s in enumerate(chips):
            a = ease((t - 0.3 - 0.13 * k) / 0.4) * out
            w = int(d.textlength(s, font=F_CHIP)) + 52
            if x + w > W - 72:
                x, yy = 72, yy + 104
            col = ORANGE if k % 3 == 0 else CYAN
            d.rounded_rectangle([x, yy + int(10 * (1 - a)), x + w, yy + 72 + int(10 * (1 - a))], 36,
                                fill=col + (int(38 * a),), outline=col + (int(220 * a),), width=2)
            text(d, (x + 26, yy + 16 + int(10 * (1 - a))), s, F_CHIP, INK, a)
            x += w + 18
    S.append((5.0, features))

    def outro(c, t, dur):
        a = ease(t / 0.8)
        lg = logo.resize((860, int(860 * logo.height / logo.width)), Image.LANCZOS)
        lg.putalpha(Image.fromarray((np.asarray(lg.getchannel("A"), float) * a).astype(np.uint8)))
        c.alpha_composite(lg, ((W - lg.width) // 2, 300))
        d = ImageDraw.Draw(c, "RGBA")
        b = ease((t - 0.6) / 0.6)
        text(d, (W // 2, 640), "github.com/sayanadhikari/picsp", font("SemiBold", 40), INK, b, anchor="ma")
        text(d, (W // 2, 705), "Open source (MIT) · 16 documented examples · verified against theory",
             font("Regular", 26), SUB, b, anchor="ma")
    S.append((4.6, outro))
    return S


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames", default=None, help="directory with <animation>_frames/ PNG folders")
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "picsp3d_showcase.mp4"))
    ap.add_argument("--preview", action="store_true", help="write one PNG per scene instead of the video")
    a = ap.parse_args()
    scenes = build_scenes(a.frames)
    total = sum(dur for dur, _ in scenes)
    if a.preview:
        t0 = 0.0
        for k, (dur, draw) in enumerate(scenes):
            c = background(t0 + dur * 0.6).convert("RGBA")
            draw(c, dur * 0.6, dur)
            progress(c, t0 + dur * 0.6, total)
            c.convert("RGB").save(os.path.splitext(a.out)[0] + f"_scene{k:02d}.png")
            t0 += dur
        return
    import imageio_ffmpeg
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    cmd = [ff, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
           "-i", "-", "-c:v", "libx264", "-preset", "slow", "-crf", "18", "-pix_fmt", "yuv420p", "-movflags",
           "+faststart", a.out]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    t0 = 0.0
    for dur, draw in scenes:
        for n in range(int(round(dur * FPS))):
            lt = n / FPS
            c = background(t0 + lt).convert("RGBA")
            draw(c, lt, dur)
            progress(c, t0 + lt, total)
            proc.stdin.write(c.convert("RGB").tobytes())
        t0 += dur
    proc.stdin.close()
    proc.wait()
    print(f"wrote {a.out} ({total:.1f} s, {W}x{H}, {FPS} fps)")


if __name__ == "__main__":
    main()
