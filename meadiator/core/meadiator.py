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
from os import getcwd
from glob import glob
from time import strftime
from collections import defaultdict
import re
import bz2
import pickle
import zipfile
from meadiator.io.meta import parse_meta, make_file_dict

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
        files_pck = "%s_files.bz2.pck" % (basename(self.base_dir))
        if exists(files_pck):
            self.files = pickle.load(bz2.BZ2File(files_pck, "r"))
            print("found existing file dictionary in %s" % (pjoin(getcwd(), files_pck)))
            for key, key_dir in [("plate", self.plate_dir)] + self.object_tups:
                print(
                    "%s found %i %s files in %s"
                    % (strftime("%H:%M:%S"), len(self.files[key]), key, key_dir)
                )
        else:
            self.files = {}
            print(
                "%s generating file lists from %s"
                % (strftime("%H:%M:%S"), self.base_dir)
            )
            self.files["plate"] = [
                x
                for x in glob(pjoin(self.plate_dir, "*", "*"))
                if basename(dirname(x)).isdigit()
                and x.split(".")[-1] in ("zip", "info")
            ]
            print(
                "%s found %i plate info files in %s"
                % (strftime("%H:%M:%S"), len(self.files["plate"]), self.plate_dir)
            )
            for key, key_dir in self.object_tups:
                self.files[key] = glob(pjoin(key_dir, r"**\*.zip"), recursive=True)
                print(
                    "%s found %i %s files in %s"
                    % (strftime("%H:%M:%S"), len(self.files[key]), key, key_dir)
                )
            pickle.dump(self.files, bz2.BZ2File(files_pck, "w"))
            print("wrote file dictionary to %s" % (pjoin(getcwd(), files_pck)))
        self.meadia = {k: {} for k, _ in self.object_tups}
        self.meadia["load_errors"] = []

    def load_objects(self):
        """
        Load MEAD objects from file paths.

        :return:
        """
        objects_pck = "%s_objects.bz2.pck" % (basename(self.base_dir))
        if exists(objects_pck):
            self.meadia = pickle.load(bz2.BZ2File(objects_pck, "r"))
            print(
                "found existing objects dictionary in %s"
                % (pjoin(getcwd(), objects_pck))
            )
        else:
            print("%s loading plate objects" % (strftime("%H:%M:%S")))
            self.meadia["plate"] = {}
            for plate_path in self.files["plate"]:
                plate_key = int(basename(plate_path).split(".")[0])
                self.meadia["plate"][plate_key] = Plate(plate_path)
            for key, key_dir in self.object_tups:
                print("%s loading %s objects" % (strftime("%H:%M:%S"), key))
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
                                    (file_path, "plate %s not in release" % (p))
                                )
                    except Exception as e:
                        self.meadia["load_errors"].append((file_path, str(e)))
            print(
                "%i files were not loaded due to read errors, see meadia['load_errors']"
                % (len(self.meadia["load_errors"]))
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
                                        "%s %s %s in plate %i info but not in release"
                                        % (otype, blk, okey, id)
                                    )
            if in_info_no_release > 0:
                print(
                    "%i runs/exps/anas are present in plate info files but were not included in the release, see meadia['load_errors']"
                    % (in_info_no_release)
                )
            # if len(self.load_errors) == 0:
            pickle.dump(self.meadia, bz2.BZ2File(objects_pck, "w"))
            print("wrote object dictionary to %s" % (pjoin(getcwd(), objects_pck)))

    def get_info(self, plate_id, return_dict=False):
        """
        Return dict of metadata for plate.

        :param plate_id: Integer plate_id.
        :return: Absolute path to info file as string.
        """
        zip_path = pjoin(self.plate_dir, str(plate_id), "%i.zip" % (plate_id))
        info_path = pjoin(dirname(zip_path), "%i.info" % (plate_id))
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
            print("no path found for %s" % (relative_path))
        elif len(found_path) > 1:
            print("multiple paths found for %s" % (relative_path))
            print("\n".join(["%i) %s" % (i, v) for i, v in enumerate(found_path)]))
        return found_path[0]

    def find_plates(
        self,
        filter_dict={
            "in_plate_id_list": [],
            "min_plate_id": None,
            "max_plate_id": None,
            "has_element": [],
            "has_deposition": [],
            "has_anneal_temp": [],
            "min_anneal_temp": None,
            "max_anneal_temp": None,
            "has_anneal_type": [],
            "has_date": [],
            "min_date": None,
            "max_date": None,
            "has_run_type": [],
            "has_exp_type": [],
            "has_ana_type": [],
            "has_ana_technique": [],
        },
        plate_list=None,
    ):
        """
        Find Plates that match criteria.

        :param filter_dict:
        :return: list of plate objects
        """
        if plate_list is None:
            plist = list(self.meadia["plate"].keys())
        else:
            plist = [x.plate_id for x in plate_list]
        min_p = (
            min(plist)
            if filter_dict["min_plate_id"] is None
            else filter_dict["min_plate_id"]
        )
        max_p = (
            max(plist)
            if filter_dict["max_plate_id"] is None
            else filter_dict["max_plate_id"]
        )
        plist = [x for x in plist if x >= min_p and x <= max_p]
        if len(filter_dict["in_plate_id_list"]) > 0:
            plist = [x for x in plist if x in filter_dict["in_plate_id_list"]]
        if len(filter_dict["has_element"]) > 0:
            plist = [
                x
                for x in plist
                if all(
                    [
                        e in self.meadia["plate"][x].elements
                        for e in filter_dict["has_element"]
                    ]
                )
            ]
        if len(filter_dict["has_run_type"]) > 0:
            plist = [
                x
                for x in plist
                if all(
                    [
                        t in self.meadia["plate"][x].run_dict.keys()
                        for t in filter_dict["has_run_type"]
                    ]
                )
            ]
        if len(filter_dict["has_exp_type"]) > 0:
            plist = [
                x
                for x in plist
                if all(
                    [
                        t in self.meadia["plate"][x].exp_dict.keys()
                        for t in filter_dict["has_exp_type"]
                    ]
                )
            ]
        if len(filter_dict["has_ana_type"]) > 0:
            plist = [
                x
                for x in plist
                if all(
                    [
                        t in self.meadia["plate"][x].ana_dict.keys()
                        for t in filter_dict["has_ana_type"]
                    ]
                )
            ]
        if len(filter_dict["has_ana_technique"]) > 0:
            plist_with_anas = []
            for p in plist:
                pobj = self.meadia["plate"][p]
                if "analyses" in vars(pobj).keys():
                    for a in pobj.analyses:
                        otype = a.split("/")[1]
                        okey = a.split("/")[-1]
                        if okey in self.meadia["ana"][otype].keys():
                            aobj = self.meadia["ana"][otype][okey]
                            if any(
                                [
                                    x in [y for y in aobj.analysis_names]
                                    for x in filter_dict["has_ana_technique"]
                                ]
                            ):
                                plist_with_anas.append(p)
            plist = set(plist_with_anas)
        return [self.meadia["plate"][p] for p in plist]

    def find_runs(
        self,
        filter_dict={
            "in_plate_id_list": [],
            "min_plate_id": None,
            "max_plate_id": None,
            "in_element_list": [],
            "has_date": [],
            "min_date": None,
            "max_date": None,
            "has_run_type": [],
        },
        plate_list=None,
    ):
        """
        Find Runs that match criteria.
        """

        return True

    def find_exps(
        self,
        filter_dict={
            "in_plate_id_list": [],
            "min_plate_id": None,
            "max_plate_id": None,
            "in_element_list": [],
            "has_date": [],
            "min_date": None,
            "max_date": None,
            "has_run_type": [],
        },
        plate_list=None,
    ):
        """
        Find Experiments that match criteria.
        """
        return True

    def find_anas(
        self,
        filter_dict={
            "in_plate_id_list": [],
            "min_plate_id": None,
            "max_plate_id": None,
            "in_element_list": [],
            "has_date": [],
            "min_date": None,
            "max_date": None,
            "has_run_type": [],
            "has_ana_techinque": [],
        },
        plate_list=None,
    ):
        """
        Find Analyses that match criteria.
        """
        return True


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

    def list_data(self):
        """
        List all data files contained in object.

        :return: list of filenames
        """
        return []

    def list_zip_contents(self):
        """
        List all contents conatined in zip.

        :return: list of filenames
        """
        return []

    def extract_file(self, file_list, target_path):
        """
        Extract filenames to target directory.

        :param file_list:
        :param target_path:
        :return:
        """
        # if isinstance(self) is "Experiment", get zip path from file_technique
        return []

    def get_zip_object(self):
        """
        Return zipfile object.

        :return:
        """
        return []

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
