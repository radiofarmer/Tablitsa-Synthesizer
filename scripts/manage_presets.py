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
import sys
import glob

BOOL = 1
CHAR = 1
WCHAR = 4
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
  DOUBLE: lambda x: np.frombuffer(x, np.float64)[0],
  CHAR: lambda x: x.decode('utf-8')
}

convert_to_bytes = {
  1: lambda x: np.array([x]).astype(np.int8).tobytes("C") if type(x) != str else x.encode("utf-8"),
  INT32: lambda x: np.array([x]).astype(np.int32).tobytes("C"),
  DOUBLE: lambda x: np.array([x]).astype(np.float64).tobytes("C"),
}

dtype_strings = {
  "bool": BOOL,
  "int32": INT32,
  "int": INT32,
  "double": DOUBLE,
  "float64": DOUBLE,
  "char": CHAR,
}

MAX_PRESET_NAME_LENGTH = 32

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
               "PresetID": INT32,
               "PresetName": [CHAR for _ in range(MAX_PRESET_NAME_LENGTH)],
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
               "VstParams": [DOUBLE for _ in range(535)]}
}


def version_string_to_hex(v_str):
  return b''.join([int(v).to_bytes(1, "little") for v in v_str.split(".")]).hex()


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
        itm_bytes = f.read(itm)
        data.update({k: convert_from_bytes[itm](itm_bytes)})
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
  v_hex = version_string_to_hex(version)
  v_int = int("0x" + v_hex, 0)
  # Get VST Parameter info
  param_filename = "v" + version
  with open(param_filename, "rb") as f:
    param_list = pickle.load(f)
  data = ReadPresetData(file_in, v_int)
  # Serialize the format specification as a list containing items of the form [PARAM_VALUE, NUM_BYTES]
  v_serial = Serialize(versions[v_int], data)

  # Actual actions:
  byte_index = 0  # Byte index for insertion actions
  if "print" in actions:  # Print all values
    print(data)
  if "modify" in actions:  # Change a value
    template = versions[v_int]
    data_bytes = bytes()
    for k in template.keys():
      param = template[k]
      if k == "VstParams":
        for p_name in param_list:
          if p_name == actions['modify']:
            data_bytes += convert_to_bytes[DOUBLE](cmdargs['value'])
          else:
            data_bytes += convert_to_bytes[DOUBLE](data[p_name])
      elif type(param) == list:
        byte_list = bytes()
        if k == actions['modify']:
          # Check data type
          if cmdargs['dtype'] == 'string':
            new_bytes = cmdargs['value'].encode('utf-8')
          else:
            vals = convert_from_string[cmdargs['dtype']](cmdargs['value'].split(''))
            new_bytes = b''.join([convert_to_bytes[cmdargs['dtype']](v) for v in vals])
          byte_list += new_bytes
          byte_list += bytes(len(param) - len(new_bytes))
        else:
          for nbytes, val in zip(param, data[k]):
            byte_list += convert_to_bytes[nbytes](val)
        data_bytes += byte_list
      else:
        if k == actions['modify']:
          data_bytes += convert_to_bytes[cmdargs['dtype']](cmdargs['value'])
        else:
          data_bytes += convert_to_bytes[param](data[k])
    with open(file_out, 'wb') as f:
      f.write(data_bytes)

  if "update" in actions:  # Change a preset to match a new format specification
    v_target = int("0x" + version_string_to_hex(actions["update"]), 0)
    data_bytes = bytes()
    for k_target in versions[v_target]:
      if k_target == "Version":
        data_bytes += convert_to_bytes[versions[v_target][k_target]](v_target)
      elif k_target not in versions[v_int]:
        param = versions[v_target][k_target]
        if type(param) == list:
          for itm in param:
            data_bytes += convert_to_bytes[itm](0)
        else:
          data_bytes += convert_to_bytes[param](0)
      elif k_target == "VstParams":
        # Get VST Parameter info for TARGET version
        param_filename_new = "v" + actions['update']
        with open(param_filename_new, "rb") as f:
          param_list_target = pickle.load(f)
        for param_name in param_list_target:
          data_bytes += convert_to_bytes[DOUBLE](data[param_name])
      else:
        param = versions[v_target][k_target]
        if type(param) == list:
          for itm, itm_target in zip(param, data[k_target]):
            data_bytes += convert_to_bytes[itm](itm_target)
        else:
          data_bytes += convert_to_bytes[param](data[k_target])
      # Save file
      if file_out[-5:].upper() != ".TPST":
        file_out += ".TPST"
      with open(file_out, "wb") as f:
        f.write(data_bytes)
    return
  if "insert-before" in actions:  # Insert a new parameter before a parameter specified by name
    # Get byte index
    for k in versions[v_int].keys():
      # Break once the chosen key has been found
      if k == actions["insert-before"]:
        break
      if type(versions[v_int][k]) == list:
        for nBytes in versions[v_int][k]:
          byte_index += nBytes
      else:
        byte_index += versions[v_int][k]
  elif "insert" in actions:
    byte_index = actions["insert"]
  if "insert" in actions or "insert-before" in actions:  # Insert a parameter before a specified byte index
    if cmdargs["value"] is None:
      print("Error: No value to be inserted was supplied.\nSpecify a value with: --value VALUE or -v VALUE")
    # Get the size and type of the data to be inserted
    dtype = cmdargs["dtype"].lower() if cmdargs["dtype"] is not None else "double"
    if cmdargs['length'] is None:
      data_size = dtype_strings[dtype]
      # Convert the data first to the specified type, then to raw bytes
      data_bytes = convert_from_string[dtype](cmdargs["value"]).to_bytes(data_size, "little")
    else:
      data_size = cmdargs['length']
      if cmdargs['dtype'] == 'string':
        data_bytes = bytearray(cmdargs['value'], 'utf-8')
        data_bytes += bytes(data_size - len(data_bytes))
      else:
        print("Unsupported data type")
        return

    if file_out[-5:].upper() != ".TPST":
      file_out += ".TPST"

    bytes_out = bytes()
    byte_count = 0
    for i, param in enumerate(v_serial):
      if byte_count == byte_index:
        bytes_out += data_bytes
      bytes_out += convert_to_bytes[param[1]](param[0])
      byte_count += param[1]

    with open(file_out, "wb") as f:
      f.write(bytes_out)


parser = argparse.ArgumentParser(description="Reformat presets for different versions")
parser.add_argument("preset", type=str)
parser.add_argument("version", type=str)
parser.add_argument("-b", "--batch", action="store_true")
parser.add_argument("-p", "--print", action="store_true")
parser.add_argument("-i", "--insert", type=int)
parser.add_argument("--insert-before", type=str)
parser.add_argument("--update", type=str)
parser.add_argument("--modify", type=str)
parser.add_argument("-v", "--value", type=str)
parser.add_argument("-d", "--dtype", type=str)
parser.add_argument("-l", "--length", type=int)
parser.add_argument("-o", "--output", type=str)

if __name__ == "__main__":
  args = parser.parse_args()
  actions_to_pass = dict()
  args_to_pass = dict()
  if args.update is not None:
    actions_to_pass.update({"update": args.update})
  elif args.modify is not None:
    actions_to_pass.update({"modify": args.modify})
    args_to_pass.update({
      "value": args.value,
      "dtype": args.dtype
    })
  elif args.insert is not None:
    actions_to_pass.update({"insert": args.insert})
    args_to_pass.update({
      "value": args.value,
      "dtype": args.dtype,
    })
  elif args.insert_before is not None:
    actions_to_pass.update({"insert-before": args.insert_before})
    args_to_pass.update({
      "value": args.value,
      "dtype": args.dtype,
      "length": args.length,
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
