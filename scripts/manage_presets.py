"""
Prints or adds data to existing Tablitsa presets, for data validation or making presets compatible with different
versions than the one they were saved in.

Usage:
  ./manage_presets -preset -version [--batch] [--print] [--insert] [--value] [--dtype] [--output]
"""

import argparse
import numpy as np
import pickle
import os
import glob

BOOL = 1
INT32 = 4
DOUBLE = 8


def string_to_int(s_int):
  if s_int[:2] == "0x":
    return int(s_int, 0)
  else:
    return int(s_int)


convert_from_string = {
  'bool': lambda s: bool(s) if s.lower() not in ["true", "false"] else s.lower() == "true",
  'int32': lambda s: string_to_int(s),
  'int': lambda s: string_to_int(s),
  'float': lambda s: float(s),
  'double': lambda s: float(s)
}

convert_from_bytes = {
  BOOL: lambda x: bool(x),
  INT32: lambda x: np.frombuffer(x, np.int32)[0],
  DOUBLE: lambda x: np.frombuffer(x, np.float64)[0]
}

convert_to_bytes = {
  BOOL: lambda x: np.array([x]).astype(np.int8).tobytes("C"),
  INT32: lambda x: np.array([x]).astype(np.int32).tobytes("C"),
  DOUBLE: lambda x: np.array([x]).astype(np.float64).tobytes("C")
}

dtype_strings = {
  "bool": BOOL,
  "int32": INT32,
  "int": INT32,
  "double": DOUBLE,
  "float64": DOUBLE
}
# "Version": INT32,

versions = {
  0x00000102: {"SequencerSteps": [DOUBLE for _ in range(16)],
               "Wavetable1": DOUBLE,
               "Wavetable2": DOUBLE,
               "SequencerQuantized": DOUBLE,
               "VoiceEffect1": INT32,
               "VoiceEffect2": INT32,
               "VoiceEffect3": INT32,
               "MasterEffect1": INT32,
               "MasterEffect1DelayTempoSync": BOOL,
               "MasterEffect2": INT32,
               "MasterEffect2DelayTempoSync": BOOL,
               "MasterEffect3": INT32,
               "MasterEffect3DelayTempoSync": BOOL,
               "VstParams": [DOUBLE for _ in range(537)]},
  0x00010100: {"Version": INT32,
               "SequencerSteps": [DOUBLE for _ in range(16)],
               "Wavetable1": DOUBLE,
               "Wavetable2": DOUBLE,
               "SequencerQuantized": DOUBLE,
               "VoiceEffect1": INT32,
               "VoiceEffect2": INT32,
               "VoiceEffect3": INT32,
               "MasterEffect1": INT32,
               "MasterEffect1DelayTempoSync": BOOL,
               "MasterEffect2": INT32,
               "MasterEffect2DelayTempoSync": BOOL,
               "MasterEffect3": INT32,
               "MasterEffect3DelayTempoSync": BOOL,
               "VstParams": [DOUBLE for _ in range(537)]}
}


def ReadPresetData(preset, version):
  v_format = versions[version]
  v_bytes = version.to_bytes(4, "little")
  data = dict()

  # Load parameter names in order for the given version
  param_filename = "v" + ".".join([str(b) for b in v_bytes[::-1]])
  with open(param_filename, "rb") as f:
    param_list = pickle.load(f)

  with open(preset, "rb") as f:
    for k in v_format.keys():
      itm = v_format[k]
      if k == "VstParams":
        for p_name, p_size in zip(param_list, v_format["VstParams"]):
          data.update({p_name: convert_from_bytes[p_size](f.read(p_size))})
      elif type(itm) == list:
        data.update({k: [convert_from_bytes[nb](f.read(nb)) for nb in itm]})
      else:
        data.update({k: convert_from_bytes[itm](f.read(itm))})
  return data


def Serialize(format_dict, data_dict):
  serial = []
  for k in format_dict.keys():
    if k == "VstParams":
      for ki, kd in enumerate([param_key for param_key in data_dict.keys() if "kParam" in param_key]):
        if "kParam" in kd:
          serial.append([data_dict[kd], format_dict["VstParams"][ki]])
    elif type(format_dict[k]) == list:
      for i, itm in enumerate(format_dict[k]):
        serial.append([data_dict[k][i], itm])
    else:
      serial.append([data_dict[k], format_dict[k]])
  return serial


def run(file_in, version, actions, file_out=None, **cmdargs):
  # Default version information
  v_hex = b''.join([int(v).to_bytes(1, "little") for v in version.split(".")]).hex()
  v_int = int("0x" + v_hex, 0)
  data = ReadPresetData(file_in, v_int)
  # Serialize the format specification as a list containing items of the form [PARAM_VALUE, NUM_BYTES]
  v_serial = Serialize(versions[v_int], data)
  # Actual actions:
  if "print" in actions:
    print(data)
  if "insert" in actions:
    if cmdargs["value"] is None:
      print("Error: No value to be inserted was supplied.\nSpecify a value with: --value VALUE or -v VALUE")
    # Get the size and type of the data to be inserted
    dtype = cmdargs["dtype"].lower() if cmdargs["dtype"] is not None else "double"
    data_size = dtype_strings[dtype]
    # Convert the data first to the specified type, then to raw bytes
    data_bytes = convert_from_string[dtype](cmdargs["value"]).to_bytes(data_size, "little")

    if file_out[-5:].upper() != ".TPST":
      file_out += ".TPST"

    bytes_out = bytes()
    for i, param in enumerate(v_serial):
      if i == actions["insert"]:
        bytes_out += data_bytes
      bytes_out += convert_to_bytes[param[1]](param[0])

    with open(file_out, "wb") as f:
      f.write(bytes_out)


parser = argparse.ArgumentParser(description="Reformat presets for different versions")
parser.add_argument("preset", type=str)
parser.add_argument("version", type=str)
parser.add_argument("-b", "--batch", action="store_true")
parser.add_argument("-p", "--print", action="store_true")
parser.add_argument("-i", "--insert", type=int)
parser.add_argument("-v", "--value", type=str)
parser.add_argument("-d", "--dtype", type=str)
parser.add_argument("-o", "--output", type=str)

if __name__ == "__main__":
  args = parser.parse_args()
  actions_to_pass = dict()
  args_to_pass = dict()
  if args.insert is not None:
    actions_to_pass.update({"insert": args.insert})
    args_to_pass.update({
      "value": args.value,
      "dtype": args.dtype,
    })
  if args.print:
    actions_to_pass.update({"print": True})

  if args.batch:
    presets = glob.glob(args.preset + "*.TPST")
    out_dir = args.output if args.output is not None else ""
    for path in presets:
      file = os.path.split(os.path.normpath(path))[-1]
      run(path, args.version, actions_to_pass, file_out=os.path.join(out_dir, file), **args_to_pass)
  else:
    output = args.output if args.output is not None else args.preset
    run(args.preset, args.version, actions_to_pass, file_out=output, **args_to_pass)