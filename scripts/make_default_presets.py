import glob
import os
import re

if __name__ == "__main__":
  preset_files = glob.glob("../resources/data/presets/*.TPST")
  preset_names = []

  p_name = re.compile(r"\/((?!.*\/).*(?=.TPST))")

  for file in preset_files:
    file_new = file.replace("\\", "/").replace(" ", "_").replace("-", "_")
    os.rename(file, file_new)
    preset_names.append(p_name.findall(file_new)[0])

  name_list = '#define PRESET_NAMES {"' + '", "'.join([pn.replace("_", " ") for pn in preset_names]) + '"}'
  id_list = '"' + '", "'.join(preset_names) + '"'
  id_var = 'constexpr char* PRESET_ID_LIST[]{' + id_list + '};'
  print(id_var)

  with open("../Presets.h", "w") as f:
    f.writelines([name_list, "\n\n", id_var])

  with open("preset_includes.txt", "w") as f:
    f.writelines([pn + ' RCDATA "' + pf.replace('\\', '/') + '"\n' for pn, pf in zip(preset_names, preset_files)])
