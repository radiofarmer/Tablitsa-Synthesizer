import glob
import os
import re
import manage_presets

if __name__ == "__main__":
  preset_files = glob.glob("../resources/data/presets/*.TPST")
  preset_names = []

  p_name = re.compile(r"\/((?!.*\/).*(?=.TPST))")

  # Get all preset file names, without full paths
  filenames = []
  for file in preset_files:
    file_new = file.replace("\\", "/").replace(" ", "_")
    preset_names.append(p_name.findall(file_new)[0])
    file_new = file_new.replace("-", "_")
    os.rename(file, file_new)
    filenames.append(file_new)

  with open('preset_names.txt', 'r') as f:
    display_names = [name[:-1] for name in f.readlines()]

  # Add names to binary data
  for name, fname in zip(display_names, filenames):
    action = {"modify": "PresetName"}
    args = {"value": name, "dtype": "string"}
    manage_presets.run(fname, "0.1.1.0", action, fname, **args)

  name_list = '#define PRESET_NAMES {"' + '", "'.join([pn.replace("_", " ") for pn in preset_names]) + '"}'
  id_list = '"' + '", "'.join(preset_names) + '"'
  id_var = 'constexpr char* PRESET_ID_LIST[]{' + id_list + '};'

  id_list_w = 'L"' + '", L"'.join(preset_names) + '"'
  id_var_w = 'constexpr wchar_t* PRESET_ID_LIST_W[]{' + id_list_w + '};'

  with open("../Presets.h", "w") as f:
    f.writelines([name_list, "\n\n", id_var, "\n\n", id_var_w, "\n"])
    f.write("\nconstexpr int N_PRESETS = {};\n".format(len(preset_names)))

  with open("preset_includes.txt", "w") as f:
    f.writelines([pn + ' RCDATA "' + os.path.abspath(pf).replace('\\', '/') + '"\n'
                  for pn, pf in zip(preset_names, preset_files)])
