import sys
import os
from PIL import Image
from pathlib import Path

if __name__ == "__main__":
  dst_path = sys.argv[1]
  src_paths = []
  for path in sys.argv[2:]:
    src_paths.append(path)

  imgdata = bytes((0, 0, 1, 0, len(src_paths), 0))
  imgs = [Image.open(p) for p in src_paths]
  offset = len(imgdata) + len(src_paths) * 16

  for im, sf in zip(imgs, src_paths):
    imgdata += bytes((im.width & 255, im.height & 255, 0, 0, 1, 0, 32, 0))
    img_size = os.path.getsize(sf)
    imgdata += img_size.to_bytes(4, "little")
    imgdata += offset.to_bytes(4, "little")
    offset += img_size

  for sf in src_paths:
    imgdata += Path(sf).read_bytes()

  Path(dst_path).write_bytes(imgdata)
