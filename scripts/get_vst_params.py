"""
Get a list of the plugin parameters which are accessible to the VST3 interface.
"""

import re
import os
import sys
import pickle

START = "enum EParams\n"
END = "};\n"


def ReadParams(filename="Tablitsa.h"):
  p = re.compile(r"\s+(.+?)[\s|,]")

  with open(os.path.join("../", filename), "r") as f:
    all_text = f.readlines()

  # Find the beginning and end of the parameter enum declaration
  start_line = all_text.index(START)
  end_line = all_text[start_line:].index("};\n")
  return [p.findall(l)[0] for l in all_text[start_line:end_line] if p.match(l)]


if __name__ == "__main__":
  param_names = ReadParams()
  with open(sys.argv[1], "wb") as f:
    pickle.dump(param_names, f)
