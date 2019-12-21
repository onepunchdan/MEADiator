"""
MEAD Internal Accessor with Thoughtful Ontology and Recall.

Provdes five object classes:
1) Meadiator
2) Plate
3) Run
4) Experiment
5) Analysis
"""
from os.path import join as pjoin, basename, dirname, exists
from os import getcwd, makedirs
from glob import glob
from time import strftime
from collections import defaultdict
import re
import bz2
import pickle
from zipfile import ZipFile
from dateutil.parser import parse
import numpy as np
import pandas as pd
from meadiator.io.meta import parse_meta, make_file_dict

PLATE_FILTER_TEMPLATE: dict = {
    "in_plate_id_list": [],
    "min_plate_id": None,
    "max_plate_id": None,
    "has_element": [],
    "has_deposition": [],
    "has_anneal_temp": [],
    "min_anneal_temp": None,
    "max_anneal_temp": None,
    "has_anneal_type": [],
    "min_date": None,
    "max_date": None,
    "has_run_type": [],
    "has_exp_type": [],
    "has_ana_type": [],
    "has_ana_technique": [],
}

RUN_FILTER_TEMPLATE: dict = {
    "in_plate_id_list": [],
    "min_plate_id": None,
    "max_plate_id": None,
    "in_element_list": [],
    "min_date": None,
    "max_date": None,
    "has_run_type": [],
}

EXP_FILTER_TEMPLATE: dict = {
    "in_plate_id_list": [],
    "min_plate_id": None,
    "max_plate_id": None,
    "in_element_list": [],
    "min_date": None,
    "max_date": None,
    "has_run_type": [],
}

ANA_FILTER_TEMPLATE: dict = {
    "in_plate_id_list": [],
    "min_plate_id": None,
    "max_plate_id": None,
    "in_element_list": [],
    "min_date": None,
    "max_date": None,
    "has_run_type": [],
    "has_ana_techinque": [],
}


class Meadiator:
    """
    MEAD Internal Accessor with Thoughtful Ontology and Recall.

    :return:
    """

    def __init__(self, base_dir):
        """
        Create meadiator object for accessing the MEAD.

        :param basedir: Path to MEAD root directory
        """
        self.log = []
        self.base_dir = base_dir
        self.run_dir = pjoin(self.base_dir, "run")
        self.plate_dir = pjoin(self.base_dir, "plate")
        self.exp_dir = pjoin(self.base_dir, "experiment")
        self.ana_dir = pjoin(self.base_dir, "analysis")
        self.object_dict = {"run": Run, "exp": Experiment, "ana": Analysis}
        self.object_tups = [
            ("run", self.run_dir),
            ("exp", self.exp_dir),
            ("ana", self.ana_dir),
        ]
        files_pck = f"{basename(self.base_dir)}_files.bz2.pck"
        if exists(files_pck):
            self.files = pickle.load(bz2.BZ2File(files_pck, "r"))
            print("found existing file dictionary in {pjoin(getcwd(), files_pck)}")
        else:
            self.files = {}
            self.log_entry(f"generating file lists from {self.base_dir}")
            self.files["plate"] = [
                x
                for x in glob(pjoin(self.plate_dir, "*", "*"))
                if basename(dirname(x)).isdigit()
                and x.split(".")[-1] in ("zip", "info")
            ]
            num_infos = len(self.files["plate"])
            self.log_entry(f"found {num_infos} plate info files in {self.plate_dir}")
            for key, key_dir in self.object_tups:
                self.files[key] = glob(pjoin(key_dir, r"**\*.zip"), recursive=True)
                self.log_entry(f"found {len(self.files[key])} {key} files in {key_dir}")
            pickle.dump(self.files, bz2.BZ2File(files_pck, "w"))
            self.log_entry("wrote file dictionary to {pjoin(getcwd(), files_pck)}")
        for key, key_dir in [("plate", self.plate_dir)] + self.object_tups:
            num_files = len(self.files[key])
            print(f"found {num_files} {key} files in {key_dir}")
        self.meadia = {k: {} for k, _ in self.object_tups}
        self.meadia["load_errors"] = []

    def log_entry(self, entry_string):
        """
        Create a log entry string in self.log and include timestamp
        """
        now = strftime("%H:%M:%S")
        self.log.append((now, entry_string))
        print(f"{now} {entry_string}")

    def load_objects(self):
        """
        Load MEAD objects from file paths.

        :return:
        """
        objects_pck = f"{basename(self.base_dir)}_objects.bz2.pck"
        if exists(objects_pck):
            self.meadia = pickle.load(bz2.BZ2File(objects_pck, "r"))
            print(
                f"found existing objects dictionary in {pjoin(getcwd(), objects_pck)}"
            )
        else:
            self.log_entry(f"loading plate objects")
            self.meadia["plate"] = {}
            for plate_path in self.files["plate"]:
                plate_key = int(basename(plate_path).split(".")[0])
                self.meadia["plate"][plate_key] = Plate(plate_path)
            for key, key_dir in self.object_tups:
                self.log_entry(f"loading {key} objects")
                for file_path in self.files[key]:
                    obj_type = file_path.replace(key_dir, "").strip("\\").split("\\")[0]
                    if obj_type not in self.meadia[key].keys():
                        self.meadia[key][obj_type] = {}
                    obj_key = basename(file_path)
                    try:
                        meadia_obj = self.object_dict[key](file_path)
                        self.meadia[key][obj_type][obj_key] = meadia_obj
                        pid = meadia_obj.plate_id
                        if isinstance(pid, str):
                            if "," in pid:
                                pids = [int(x.strip()) for x in pid.split(",")]
                        elif isinstance(pid, list):
                            pids = [int(x) for x in pid]
                        else:
                            pids = [pid]
                        for p in pids:
                            if p in self.meadia["plate"].keys():
                                if key == "run":
                                    if (
                                        obj_type
                                        not in self.meadia["plate"][p].run_dict.keys()
                                    ):
                                        self.meadia["plate"][p].run_dict[obj_type] = {}
                                    self.meadia["plate"][p].run_dict[obj_type][
                                        obj_key
                                    ] = meadia_obj
                                elif key == "exp":
                                    if (
                                        obj_type
                                        not in self.meadia["plate"][p].exp_dict.keys()
                                    ):
                                        self.meadia["plate"][p].exp_dict[obj_type] = {}
                                    self.meadia["plate"][p].exp_dict[obj_type][
                                        obj_key
                                    ] = meadia_obj
                                elif key == "ana":
                                    if (
                                        obj_type
                                        not in self.meadia["plate"][p].ana_dict.keys()
                                    ):
                                        self.meadia["plate"][p].ana_dict[obj_type] = {}
                                    self.meadia["plate"][p].ana_dict[obj_type][
                                        obj_key
                                    ] = meadia_obj
                            else:
                                self.meadia["load_errors"].append(
                                    (file_path, f"plate {p} not in release")
                                )
                    except Exception as e:
                        self.meadia["load_errors"].append((file_path, str(e)))
            num_errors = len(self.meadia["load_errors"])
            self.log_entry(
                f"{num_errors} files were not loaded due to read errors, see meadia['load_errors']"
            )
            in_info_no_release = 0
            for plate_path in self.files["plate"]:  # propogate meta data
                plate_meta = parse_meta(plate_path)
                id = plate_meta["plate_id"]
                elements = self.meadia["plate"][id].elements
                ann_temp = self.meadia["plate"][id].anneal_temp
                ann_type = self.meadia["plate"][id].anneal_type
                for block in ["runs", "experiments", "analyses"]:
                    blk = block[:3]
                    if block in plate_meta.keys():
                        # update date, elements, anneal_temp, anneal_type
                        if isinstance(plate_meta[block], dict):
                            for k, blkd in plate_meta[block].items():
                                otype = blkd["path"].split("/")[1]
                                okey = blkd["path"].split("/")[-1]
                                if okey in self.meadia[blk][otype].keys():
                                    self.meadia[blk][otype][okey].elements = elements
                                    self.meadia[blk][otype][okey].anneal_temp = ann_temp
                                    self.meadia[blk][otype][okey].anneal_type = ann_type
                                    self.meadia[blk][otype][okey].anneal_temp = ann_temp
                                    self.meadia[blk][otype][okey].date = blkd[
                                        "created_at"
                                    ]
                                    if blk == "run":
                                        # update machine, file_count
                                        if "machine" in blkd.keys():
                                            self.meadia[blk][otype][
                                                okey
                                            ].machine = blkd["machine"]
                                        if "description" in blkd.keys():
                                            self.meadia[blk][otype][okey].file_count = (
                                                blkd["description"]
                                                .split("containing ")[1]
                                                .split(" files")[0]
                                            )
                                    elif blk == "exp":
                                        # update runs
                                        if "run_paths" in blkd.keys():
                                            self.meadia[blk][otype][okey].runs = [
                                                self.meadia["runs"][otype][basename(p)]
                                                for p in blkd["run_paths"]
                                            ]
                                    elif blk == "ana":
                                        # update experiments
                                        if "experiment_path" in blkd.keys():
                                            self.meadia[blk][otype][
                                                okey
                                            ].experiment = self.meadia["experiments"][
                                                otype
                                            ][
                                                basename(blkd["experiment_path"])
                                            ]
                                else:
                                    in_info_no_release += 1
                                    self.meadia["load_errors"].append(
                                        f"{otype} {blk} {okey} in plate {id} info but not in release"
                                    )
            if in_info_no_release > 0:
                self.log_entry(
                    f"{in_info_no_release} runs/exps/anas are present in plate info files but were not included in the release, see meadia['load_errors']"
                )
            # if len(self.load_errors) == 0:
            pickle.dump(self.meadia, bz2.BZ2File(objects_pck, "w"))
            self.log_entry(f"wrote object dictionary to {pjoin(getcwd(), objects_pck)}")

    def get_info(self, plate_id, return_dict=False):
        """
        Return dict of metadata for plate.

        :param plate_id: Integer plate_id.
        :return: Absolute path to info file as string.
        """
        zip_path = pjoin(self.plate_dir, str(plate_id), f"{plate_id}.zip")
        info_path = pjoin(dirname(zip_path), f"{plate_id}.info")
        if exists(zip_path):
            file_path = zip_path
        elif exists(info_path):
            file_path = info_path
        else:
            file_path = ""
        return parse_meta(file_path) if return_dict else file_path

    def get_path(self, dir_key, relative_path):
        """
        Return absolute path from relative path and directory.

        :param dir_key:
        :param relative_path:
        :return: Absolute file path as string.
        """
        dir_dict = {"run": self.run_dir, "exp": self.exp_dir, "ana": self.ana_dir}
        found_path = glob(
            pjoin(
                dir_dict[dir_key],
                ".".join(relative_path.split(".")[:2]).lstrip("/") + ".*",
            )
        )
        if len(found_path) == 0:
            found_path = [""]
            self.log_entry(f"no path found for {relative_path}")
        elif len(found_path) > 1:
            self.log_entry(f"multiple paths found for {relative_path}")
            for i, v in enumerate(found_path):
                self.log_entry(f"{i}) {v}")
        return found_path[0]

    def find_plates(self, filter_dict=PLATE_FILTER_TEMPLATE, plate_list=None):
        """
        Find Plates that match criteria.

        :param filter_dict:
        :return: list of plate objects
        """
        if plate_list is None:
            plist = list(self.meadia["plate"].keys())
        else:
            plist = [x.plate_id for x in plate_list]

        # filter by plate id
        if "min_plate_id" in filter_dict.keys():
            min_p = filter_dict["min_plate_id"]
        else:
            min_p = min(plist)

        if "max_plate_id" in filter_dict.keys():
            max_p = filter_dict["max_plate_id"]
        else:
            max_p = max(plist)
        plist = [x for x in plist if min_p <= x <= max_p]

        if "in_plate_id_list" in filter_dict.keys():
            pids = filter_dict["in_plate_id_list"]
            if isinstance(pids, int):
                pids = [pids]
            plist = [x for x in plist if x in pids]

        if "has_element" in filter_dict.keys():
            els = filter_dict["has_element"]
            if isinstance(els, str):
                els = [els]
            plist = [x for x in plist if all([e in self.meadia["plate"][x].elements for e in els])]

        if "has_run_type" in filter_dict.keys():
            rtypes = filter_dict["has_run_type"]
            if isinstance(rtypes, str):
                rtypes = [rtypes]
            plist = [x for x in plist if all([t in self.meadia["plate"][x].run_dict.keys() for t in rtypes])]

        if "has_exp_type" in filter_dict.keys():
            etypes = filter_dict["has_exp_type"]
            if isinstance(etypes, str):
                etypes = [etypes]
            plist = [x for x in plist if all([t in self.meadia["plate"][x].exp_dict.keys() for t in etypes])]

        if "has_ana_type" in filter_dict.keys():
            atypes = filter_dict["has_ana_type"]
            if isinstance(atypes, str):
                atypes = [atypes]
            plist = [x for x in plist if all([t in self.meadia["plate"][x].ana_dict.keys() for t in atypes])]

        if "has_ana_technique" in filter_dict.keys():
            atechs = filter_dict["has_ana_technique"]
            if isinstance(atechs, str):
                atechs = [atechs]
            plist_with_anas = []
            for p in plist:
                pobj = self.meadia["plate"][p]
                if "analyses" in vars(pobj).keys():
                    for a in pobj.analyses:
                        otype = a.split("/")[1]
                        okey = a.split("/")[-1]
                        if okey in self.meadia["ana"][otype].keys():
                            aobj = self.meadia["ana"][otype][okey]
                            if any([all([x in y for x in atechs]) for y in aobj.analysis_names]):
                                plist_with_anas.append(p)
            plist = set(plist_with_anas)

        if "min_date" in filter_dict.keys():
            min_ts = filter_dict['min_date']
            if isinstance(min_ts, str):
                min_ts = parse(min_ts)
            plist = [p for p in plist if parse(self.meadia["plate"][p].date)>=min_ts]

        if "max_date" in filter_dict.keys():
            max_ts = filter_dict['max_date']
            if isinstance(max_ts, str):
                max_ts = parse(max_ts)
            plist = [p for p in plist if parse(self.meadia["plate"][p].date)<=max_ts]

        return [self.meadia["plate"][p] for p in plist]

    def find_runs(self, filter_dict=RUN_FILTER_TEMPLATE, plate_list=None):
        """
        Find Runs that match criteria.
        """
        if plate_list is None:
            plist = list(self.meadia["plate"].keys())
        else:
            plist = [x.plate_id for x in plate_list]
        # filter by plate id
        if "min_plate_id" in filter_dict.keys():
            min_p = filter_dict["min_plate_id"]
        else:
            min_p = min(plist)

        if "max_plate_id" in filter_dict.keys():
            max_p = filter_dict["max_plate_id"]
        else:
            max_p = max(plist)
        plist = [x for x in plist if min_p <= x <= max_p]

        if "in_plate_id_list" in filter_dict.keys():
            pids = filter_dict["in_plate_id_list"]
            if isinstance(pids, int):
                pids = [pids]
            plist = [x for x in plist if x in pids]

        if "has_element" in filter_dict.keys():
            els = filter_dict["has_element"]
            if isinstance(els, str):
                els = [els]
            plist = [x for x in plist if all([e in self.meadia["plate"][x].elements for e in els])]

        rlist = []
        if "has_run_type" in filter_dict.keys():
            rtype = filter_dict["has_run_type"]
            if isinstance(rtype, str):
                rtype = [rtype]
        else:
            rtype = list(self.meadia["run"].keys())
        for t in rtype:
            if t not in self.meadia["run"].keys():
                print("Run type {t} does not exist in MEAD data set")
                continue
            for _, run in self.meadia["run"][t].items():
                pid = run.plate_id
                if isinstance(pid, int):
                    pid = [pid]
                if any([p in plist for p in pid]):
                    rlist.append(run)

        if "min_date" in filter_dict.keys():
            min_ts = filter_dict['min_date']
            if isinstance(min_ts, str):
                min_ts = parse(min_ts)
            rlist = [r for r in rlist if parse(r.date)>=min_ts]

        if "max_date" in filter_dict.keys():
            max_ts = filter_dict['max_date']
            if isinstance(max_ts, str):
                max_ts = parse(max_ts)
            rlist = [r for r in rlist if parse(r.date)<=max_ts]

        return rlist

    def find_exps(self, filter_dict=EXP_FILTER_TEMPLATE, plate_list=None):
        """
        Find Experiments that match criteria.
        """
        if plate_list is None:
            plist = list(self.meadia["plate"].keys())
        else:
            plist = [x.plate_id for x in plate_list]
        # filter by plate id
        if "min_plate_id" in filter_dict.keys():
            min_p = filter_dict["min_plate_id"]
        else:
            min_p = min(plist)

        if "max_plate_id" in filter_dict.keys():
            max_p = filter_dict["max_plate_id"]
        else:
            max_p = max(plist)
        plist = [x for x in plist if min_p <= x <= max_p]

        if "in_plate_id_list" in filter_dict.keys():
            pids = filter_dict["in_plate_id_list"]
            if isinstance(pids, int):
                pids = [pids]
            plist = [x for x in plist if x in pids]

        if "has_element" in filter_dict.keys():
            els = filter_dict["has_element"]
            if isinstance(els, str):
                els = [els]
            plist = [x for x in plist if all([e in self.meadia["plate"][x].elements for e in els])]

        elist = []
        if "has_run_type" in filter_dict.keys():
            etype = filter_dict["has_run_type"]
            if isinstance(etype, str):
                etype = [etype]
        else:
            etype = list(self.meadia["exp"].keys())
        for t in etype:
            if t not in self.meadia["exp"].keys():
                print("Experiment type {t} does not exist in MEAD data set")
                continue
            for _, exp in self.meadia["exp"][t].items():
                pid = exp.plate_id
                if isinstance(pid, int):
                    pid = [pid]
                if any([p in plist for p in pid]):
                    elist.append(exp)

        if "min_date" in filter_dict.keys():
            min_ts = filter_dict['min_date']
            if isinstance(min_ts, str):
                min_ts = parse(min_ts)
            elist = [e for e in elist if parse(e.date)>=min_ts]

        if "max_date" in filter_dict.keys():
            max_ts = filter_dict['max_date']
            if isinstance(max_ts, str):
                max_ts = parse(max_ts)
            elist = [e for e in elist if parse(e.date)<=max_ts]

        return elist

    def find_anas(self, filter_dict=ANA_FILTER_TEMPLATE, plate_list=None):
        """
        Find Analyses that match criteria.
        """
        if plate_list is None:
            plist = list(self.meadia["plate"].keys())
        else:
            plist = [x.plate_id for x in plate_list]
        # filter by plate id
        if "min_plate_id" in filter_dict.keys():
            min_p = filter_dict["min_plate_id"]
        else:
            min_p = min(plist)

        if "max_plate_id" in filter_dict.keys():
            max_p = filter_dict["max_plate_id"]
        else:
            max_p = max(plist)
        plist = [x for x in plist if min_p <= x <= max_p]

        if "in_plate_id_list" in filter_dict.keys():
            pids = filter_dict["in_plate_id_list"]
            if isinstance(pids, int):
                pids = [pids]
            plist = [x for x in plist if x in pids]

        if "has_element" in filter_dict.keys():
            els = filter_dict["has_element"]
            if isinstance(els, str):
                els = [els]
            plist = [x for x in plist if all([e in self.meadia["plate"][x].elements for e in els])]

        alist = []
        if "has_ana_type" in filter_dict.keys():
            atype = filter_dict["has_ana_type"]
            if isinstance(atype, str):
                atype = [atype]
        else:
            atype = list(self.meadia["ana"].keys())
        for t in atype:
            if t not in self.meadia["ana"].keys():
                print("Analysis type {t} does not exist in MEAD data set")
                continue
            for _, ana in self.meadia["ana"][t].items():
                pid = ana.plate_id
                if isinstance(pid, int):
                    pid = [pid]
                if any([p in plist for p in pid]):
                    alist.append(ana)

        if "min_date" in filter_dict.keys():
            min_ts = filter_dict['min_date']
            if isinstance(min_ts, str):
                min_ts = parse(min_ts)
            alist = [a for a in alist if parse(a.date)>=min_ts]

        if "max_date" in filter_dict.keys():
            max_ts = filter_dict['max_date']
            if isinstance(max_ts, str):
                max_ts = parse(max_ts)
            alist = [a for a in alist if parse(a.date)<=max_ts]

        return alist


class Plate:
    """
    MEAD Plate object.

    fill in description
    """

    def __init__(self, info_path):
        """
        Plate object returned by Meadiator.

        :param info_path:
        :return:
        """
        tmpd = defaultdict(str, parse_meta(info_path))
        self.path = info_path
        self.plate_id = tmpd["plate_id"]
        self.date = tmpd["created_at"]
        self.lineage = tmpd["lineage"]
        desc = tmpd["description"]
        desc_elstring = desc.replace("Material library with ", "").split()[0]
        desc_elstring = re.sub("([A-Za-z])([A-Z])", "\\1,\\2", desc_elstring)
        desc_els = desc_elstring.split(",")
        if isinstance(desc_els, str):
            desc_els = [desc_els]
        self.elements = list(set(desc_els))
        self.deposition_method = desc.split("deposited by ")[-1].split()[0]
        self.substrate = tmpd["substrate"]
        self.map = tmpd["screening_map_id"]
        if "annealed at" in desc:
            self.anneal_temp = float(
                desc.split("annealed at ")[-1].split()[0].replace("C", "").strip()
            )
        else:
            self.anneal_temp = 0
        if "to add" in desc:
            self.anneal_type = desc.split("to add ")[-1].split()[0]
        else:
            self.anneal_type = ""
        if " on " in desc:
            self.anneal_date = desc.split(" on ")[-1].strip()
        else:
            self.anneal_date = ""
        self.run_dict = {}
        self.exp_dict = {}
        self.ana_dict = {}
        self.runs = []
        if "runs" in tmpd.keys():
            if isinstance(tmpd["runs"], dict):
                for rund in tmpd["runs"].values():
                    self.runs.append(rund["path"])
        self.experiments = []
        if "experiments" in tmpd.keys():
            if isinstance(tmpd["experiments"], dict):
                for expd in tmpd["experiments"].values():
                    self.experiments.append(expd["path"])
        self.analyses = []
        if "analyses" in tmpd.keys():
            if isinstance(tmpd["analyses"], dict):
                for anad in tmpd["analyses"].values():
                    self.analyses.append(anad["path"])
        del tmpd


class Meadia:
    """Meadia parent class common methods."""

    def __init__(self):
        self.files = {}
        self.path = ""

    def list_data(self):
        """
        List all data files contained in object.

        :return: list of filenames
        """
        data_list = []
        for _, techd in self.files.items():
            for ftype, filed in techd.items():
                if "_files" in ftype:
                    data_list += list(filed.keys())
        return data_list

    def list_zip_contents(self):
        """
        List all contents conatined in zip.

        :return: list of filenames
        """
        with ZipFile(self.path) as z:
            file_list = z.namelist()
        return file_list

    def extract_file(self, file_list, target_path):
        """
        Extract filenames to target directory.

        :param file_list:
        :param target_path:
        :return:
        """
        # if isinstance(self) is "Experiment", get zip path from file_technique
        zfile_dict = []
        for _, techd in self.files.items():
            for ftype, filed in techd.items():
                if "_files" in ftype:
                    for k in filed.keys():
                        zfile_dict[k] = filed[k]["source_zip"]

        for f in file_list:
            try:
                makedirs(target_path, exist_ok=True)
                print("")
            except Exception:
                print("")

    def read_lines(self, data_file):
        """
        Basic read_lines from file. Convenience wrapper for zipped files.

        :param data_file:
        :return:
        """
        return []

    def read_nparray(self, data_file, skip):
        """
        Return numpy array of data. Intended for image files
        """
        # if isinstance(self) is "Run", load from first numeric to last numeric line
        # check if data_file is image
        return []

    def read_dataframe(self, data_file, skip, cols):
        """
        """
        # if isinstance(self) is "Run", load from first numeric to last numeric line
        # check if data_file is image -> can't load
        return []


class Run(Meadia):
    """
    MEAD Run object.

    fill in description
    """

    def __init__(self, run_path):
        """
        Run object returned by Meadiator.

        :param run_path:
        :return:
        """
        tmpd = defaultdict(str, parse_meta(run_path))
        self.path = tmpd["file_path"]
        self.date = ""
        self.type = tmpd["experiment_type"]
        self.plate_id = tmpd["plate_id"]
        self.machine = ""
        self.elements = []
        self.anneal_temp = 0
        self.anneal_type = ""
        self.file_count = 0
        common_keys = ["plate_id", "experiment_type"]
        file_keys = [k for k in tmpd.keys() if "files_technique__" in k]
        self.files = {}
        for key in file_keys:
            file_tech = key.split("__")[-1]
            self.files[file_tech] = {}
            for file_group in tmpd[key].keys():
                if isinstance(tmpd[key][file_group], dict):
                    file_dict = make_file_dict(tmpd[key][file_group])
                    for v in file_dict.values():
                        v.update(source_zip=self.path)
                    self.files[file_tech][file_group] = file_dict
                else:
                    continue
        param_keys = [k for k in tmpd.keys() if "params__" in k or k == "parameters"]
        self.techs = list(self.files.keys())
        self.tech_params = {}
        tech_param_list = []
        for key in param_keys:
            self.tech_params[key.split("__")[-1]] = tmpd[key]
            tech_param_list += list(tmpd[key].keys())
        self.tech_param_groups = list(self.tech_params.keys())
        self.tech_param_keys = list(set(tech_param_list))
        other_keys = [
            k for k in tmpd.keys() if k not in common_keys + file_keys + param_keys
        ]
        self.root_params = {k: tmpd[k] for k in other_keys}
        self.root_keys = other_keys
        del tmpd


class Experiment(Meadia):
    """
    MEAD Experiment object.

    Fill in description
    """

    def __init__(self, exp_path):
        """
        Experiment object returned by Meadiator.

        :param exp_path:
        :return:
        """
        tmpd = defaultdict(str, parse_meta(exp_path))
        self.path = tmpd["file_path"]
        self.date = ""
        self.type = tmpd["experiment_type"]
        self.plate_id = tmpd["plate_ids"]
        if isinstance(self.plate_id, str):
            self.plate_id = [int(x.strip()) for x in self.plate_id.split(",")]
        self.elements = []
        self.anneal_temp = 0
        self.anneal_type = ""
        run_keys = [k for k in tmpd.keys() if k.startswith("run__")]
        self.run_count = len(run_keys)
        self.run_paths = [tmpd[k]["run_path"] for k in run_keys]
        self.files = {}
        for run_key in run_keys:
            file_keys = [k for k in tmpd[run_key].keys() if "files_technique__" in k]
            for key in file_keys:
                file_tech = key.split("__")[-1]
                if file_tech not in self.files.keys():
                    self.files[file_tech] = {}
                for file_group in tmpd[run_key][key].keys():
                    if file_group not in self.files[file_tech].keys():
                        self.files[file_tech][file_group] = {}
                    if isinstance(tmpd[run_key][key][file_group], dict):
                        file_dict = make_file_dict(tmpd[run_key][key][file_group])
                        for v in file_dict.values():
                            v.update(source_zip=self.path)
                        self.files[file_tech][file_group].update(file_dict)
                    else:
                        continue
        self.techs = list(self.files.keys())
        del tmpd
        self.runs = []

    def get_runs(self, meadia_dict):
        """
        Populate runs attribute with MEAD Runs.

        :param meadia_dict:
        :return:
        """
        for p in self.run_paths:
            run_type = p.strip("/").split("/")[0]
            filename = p.split("/")[-1]
            self.runs.append(meadia_dict["run"][run_type][filename])


class Analysis(Meadia):
    """
    MEAD Analysis object.

    Fill in description
    """

    def __init__(self, ana_path):
        """
        Experiment object returned by Meadiator.

        :param exp_path:
        :return:
        """
        tmpd = defaultdict(str, parse_meta(ana_path))
        self.path = tmpd["file_path"]
        self.date = ""
        self.type = tmpd["experiment_type"]
        self.plate_id = tmpd["plate_ids"]
        if isinstance(self.plate_id, str):
            self.plate_id = [int(x.strip()) for x in self.plate_id.split(",")]
        self.elements = []
        self.anneal_temp = 0
        self.anneal_type = ""
        self.experiment_path = tmpd["experiment_path"]
        self.analysis_params = {}
        self.files = {}
        ana_keys = [k for k in tmpd.keys() if k.startswith("ana__")]
        self.ana_count = len(ana_keys)
        for ana_key in ana_keys:
            self.analysis_params[ana_key] = {"name": tmpd[ana_key]["name"]}
            if "parameters" in tmpd[ana_key].keys():
                self.analysis_params[ana_key].update(tmpd[ana_key]["parameters"])
            file_keys = [k for k in tmpd[ana_key].keys() if "files_" in k]
            if "technique" in tmpd[ana_key].keys():
                file_tech = tmpd[ana_key]["technique"]
            elif "analysis_general_type" in tmpd[ana_key].keys():
                if tmpd[ana_key]["analysis_general_type"] == "process_fom":
                    sourced = tmpd[tmpd[ana_key]["parameters"]["select_ana"]]
                    while "technique" not in sourced.keys():
                        sourced = tmpd[sourced["parameters"]["select_ana"]]
                    file_tech = sourced["technique"]
                else:
                    file_tech = "no_technique"
            else:
                file_tech = "no_technique"
            for key in file_keys:
                if file_tech not in self.files.keys():
                    self.files[file_tech] = {}
                for file_group in tmpd[ana_key][key].keys():
                    if file_group not in self.files[file_tech].keys():
                        self.files[file_tech][file_group] = {}
                    if isinstance(tmpd[ana_key][key][file_group], dict):
                        file_dict = make_file_dict(tmpd[ana_key][key][file_group])
                        for v in file_dict.values():
                            v.update(source_zip=self.path)
                        self.files[file_tech][file_group].update(file_dict)
                    else:
                        continue
        self.techs = list(self.files.keys())
        self.analysis_names = [self.analysis_params[k]["name"] for k in ana_keys]
        del tmpd
        self.experiment = None

    def get_experiment(self, meadia_dict):
        """
        Populate experiment_path attribute with MEAD Experiment.

        :param meadia_dict:
        :return:
        """
        exp_type = self.experiment_path.strip("/").split("/")[0]
        filename = self.experiment_path.split("/")[-1]
        self.experiment = meadia_dict["exp"][exp_type][filename]
