"""GIF 保存の共通処理: 全コマで同じパレット（色の表）を使い、コマごとの色のちらつきを防ぐ"""
from PIL import Image


def save_gif(frames_rgb, path, duration=70, colors=192):
    W, H = frames_rgb[0].size
    pick = frames_rgb[::max(1, len(frames_rgb) // 6)]
    mosaic = Image.new("RGB", (W, H * len(pick)))
    for k, f in enumerate(pick):
        mosaic.paste(f.convert("RGB"), (0, H * k))
    pal = mosaic.quantize(colors=colors, method=Image.Quantize.MEDIANCUT)
    out = [f.convert("RGB").quantize(palette=pal, dither=Image.Dither.NONE) for f in frames_rgb]
    out[0].save(path, save_all=True, append_images=out[1:], duration=duration, loop=0, optimize=False, disposal=1)
