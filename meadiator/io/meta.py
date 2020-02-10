from os.path import dirname, basename
from collections import defaultdict

import zipfile

def isfloat(value):
    try:
        float(value)
        return True
    except Exception:
        return False

def parse_meta(file_path):
    """
    Read metadata from file.

    Compatible filetypes are ".rcp" ".exp" ".ana" ".info" ".zip"
    :param file_path: Full path to metadata file.
    :return: Metadata parsed into dict (often nested).
    """

    def tab_level(any_string):
        """
        Count number of leading tabs in a string.

        :param any_string:
        :return:
        """
        return (len(any_string) - len(any_string.lstrip("    "))) / 4

    tup_list = []

    try:
        if file_path.endswith(".zip"):
            if "analysis" in dirname(file_path):
                ext = ".ana"
            elif "experiment" in dirname(file_path):
                ext = ".exp"
            elif "run" in dirname(file_path):
                ext = ".rcp"
            elif "plate" in dirname(file_path):
                ext = ".info"
            meta_file = basename(file_path).split(".copied")[0]
            if ext not in [".ana", ".exp", ".rcp"]:
                meta_file = meta_file.split("-")[0].split(".zip")[0]
            meta_file += ext
            archive = zipfile.ZipFile(file_path, "r")
            with archive.open(meta_file, "r") as f:
                for l in f:
                    if l.decode("ascii").strip() != "":
                        k, v = l.decode("ascii").split(":", 1)
                        lvl = tab_level(l.decode("ascii"))
                        tup_list.append((lvl, k.strip(), v.strip()))
        else:
            with open(file_path, "r") as f:
                for l in f:
                    if l.strip() != "":
                        k, v = l.split(":", 1)
                        lvl = tab_level(l)
                        tup_list.append((lvl, k.strip(), v.strip()))
    except:
        print(f"Could not read metafile in {file_path}")
        return {"file_path": file_path}

    def build_dict(tups, clvl=0):
        sub_dict = {}
        for j, tup in enumerate(tups):
            lvl, key, val = tup
            lvl -= clvl
            if lvl != 0:
                continue
            if val == "":
                root_lvls = [
                    k for k, tup in enumerate(tups) if tup[0] - clvl == 0 and k > j
                ]
                next_rootlvl = min(root_lvls) if len(root_lvls) > 0 else len(tups)
                sub_dict[key] = build_dict(tups[j + 1 : next_rootlvl], clvl + 1)
            else:
                if isfloat(val):
                    val = float(val)
                    if val.is_integer():
                        val = int(val)
                sub_dict[key] = val
        return sub_dict

    final_dict = build_dict(tup_list)

    return final_dict



def make_file_dict(d):
    filed = {}
    for k, v in d.items():
        fn = k
        vlist = v.strip(";").split(";")
        if len(vlist) == 5:
            metad = {
                mk: mv
                for mk, mv in zip(
                    ("file_type", "names", "skip", "stopind", "sample_no"),
                    vlist
                )
            }
        elif len(vlist) == 4:
            metad = {
                mk: mv
                for mk, mv in zip(
                    ("file_type", "names", "skip", "stopind"),
                    vlist
                )
            }
        elif len(vlist) == 2:
            if vlist[-1].isdigit():
                metad = {
                    mk: mv
                    for mk, mv in zip(
                        ("file_type", "sample_no"),
                        vlist
                    )
                }
            else:
                metad = {
                    mk: mv
                    for mk, mv in zip(
                        ("file_type", "names"),
                        vlist
                    )
                }
        else:
            metad = {"file_type": vlist[0]}
        filed[fn] = metad
    return filed
