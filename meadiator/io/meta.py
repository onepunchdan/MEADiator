from os.path import dirname, basename

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
    dict_list = []

    def tab_level(any_string):
        """
        Count number of leading tabs in a string.

        :param any_string:
        :return:
        """
        return (len(any_string) - len(any_string.lstrip("    "))) / 4

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
                        dict_list.append(
                            {"name": k.strip(), "value": v.strip(), "level": lvl}
                        )
        else:
            with open(file_path, "r") as f:
                for l in f:
                    if l.strip() != "":
                        k, v = l.split(":", 1)
                        lvl = tab_level(l)
                        dict_list.append(
                            {"name": k.strip(), "value": v.strip(), "level": lvl}
                        )
    except:
        print("Could not read metafile in %s" % (file_path))
        return {'file_path': file_path}

    def ttree_to_json(ttree, level=0):
        result = {}
        for i in range(0, len(ttree)):
            cn = ttree[i]
            try:
                nn = ttree[i + 1]
            except:
                nn = {"level": -1}
            if cn["level"] > level:
                continue
            if cn["level"] < level:
                return result
            if nn["level"] == level:
                dict_insert_or_append(result, cn["name"], cn["value"])
            elif nn["level"] > level:
                rr = ttree_to_json(ttree[i + 1 :], level=nn["level"])
                dict_insert_or_append(result, cn["name"], rr)
            else:
                dict_insert_or_append(result, cn["name"], cn["value"])
                return result
        return result

    def dict_insert_or_append(a_dict, key, val):
        """
        Create or append to existing key.

        Insert a value in dict at key if one does not exist. Otherwise, convert
        value to list and append.
        :param a_dict:
        :param key:
        :param val:
        :return:
        """
        if isfloat(val):
            val = float(val)
            if val.is_integer():
                val = int(val)
        if key in a_dict:
            if type(a_dict[key]) != list:
                a_dict[key] = [a_dict[key]]
            a_dict[key].append(val)
        else:
            a_dict[key] = val

    parsed = ttree_to_json(dict_list)
    parsed.update({"file_path": file_path})

    return parsed


def make_file_dict(d):
    filed = {}
    for k, v in d.items():
        fn = k
        vlist = v.strip(";").split(";")
        if len(vlist) == 5:
            metad = {
                mk: mv
                for mk, mv in zip(
                    ("file_type", "names", "skip", "data_len", "sample_no"),
                    vlist
                )
            }
        elif len(vlist) == 4:
            metad = {
                mk: mv
                for mk, mv in zip(
                    ("file_type", "names", "skip", "data_len"),
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
