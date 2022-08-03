Search.setIndex({"docnames": ["HPC_submissions", "TIES_MD", "TIES_MD.eng_scripts", "TIES_MD.eng_scripts.cfg_scripts", "TIES_MD.eng_scripts.namd_sub", "TIES_MD.eng_scripts.namd_sub_split", "TIES_MD.eng_scripts.openmm_sub_split", "TIES_MD.openmmtools", "TIES_MD.ties_analysis", "TIES_MD.ties_analysis.engines", "TIES_MD.ties_analysis.methods", "binding_free_energies", "index", "install", "modules", "parallelization", "theory", "tutorial"], "filenames": ["HPC_submissions.rst", "TIES_MD.rst", "TIES_MD.eng_scripts.rst", "TIES_MD.eng_scripts.cfg_scripts.rst", "TIES_MD.eng_scripts.namd_sub.rst", "TIES_MD.eng_scripts.namd_sub_split.rst", "TIES_MD.eng_scripts.openmm_sub_split.rst", "TIES_MD.openmmtools.rst", "TIES_MD.ties_analysis.rst", "TIES_MD.ties_analysis.engines.rst", "TIES_MD.ties_analysis.methods.rst", "binding_free_energies.rst", "index.rst", "install.rst", "modules.rst", "parallelization.rst", "theory.rst", "tutorial.rst"], "titles": ["HPC Submission scripts", "TIES_MD package", "TIES_MD.eng_scripts package", "TIES_MD.eng_scripts.cfg_scripts package", "TIES_MD.eng_scripts.namd_sub package", "TIES_MD.eng_scripts.namd_sub_split package", "TIES_MD.eng_scripts.openmm_sub_split package", "TIES_MD.openmmtools package", "TIES_MD.ties_analysis package", "TIES_MD.ties_analysis.engines package", "TIES_MD.ties_analysis.methods package", "Binding Free Energy Tutorial", "Welcome to the TIES MD documentation.", "Installation", "TIES_MD", "Parallelization", "Theory", "Tutorial"], "terms": {"here": [0, 1, 13, 15, 16, 17], "we": [0, 1, 8, 9, 10, 11, 15, 16, 17], "provid": [0, 1, 13, 15, 16, 17], "some": [0, 1, 11, 13, 15, 16, 17], "exampl": [0, 1, 11, 13, 15, 16, 17], "variou": 0, "system": [0, 1, 7, 11, 12, 13, 15, 16, 17], "ti": [0, 8, 9, 11, 14, 16, 17], "md": [0, 11, 15, 16, 17], "attempt": [0, 13], "automat": [0, 17], "write": [0, 1, 11, 13, 15, 17], "sensibl": 0, "namd2": [0, 1, 15, 17], "target": [0, 1, 13, 16, 17], "archer": [0, 15], "2": [0, 9, 11, 13, 15, 16, 17], "summit": [0, 1, 11, 15], "In": [0, 11, 13, 15, 17], "gener": [0, 8, 9, 12, 15, 16, 17], "user": [0, 1, 11, 13, 16, 17], "can": [0, 1, 7, 10, 11, 12, 13, 15, 16, 17], "make": [0, 1, 10, 11, 15], "own": [0, 11, 17], "whichev": 0, "cluster": [0, 11], "thei": [0, 1, 11], "prefer": 0, "To": [0, 11, 13, 15, 16, 17], "aid": 0, "expos": 0, "option": [0, 1, 7, 8, 9, 11, 15, 16, 17], "sub_head": [0, 11], "sub_run_lin": [0, 11], "The": [0, 1, 7, 11, 13, 15, 16, 17], "string": [0, 1, 8, 9, 10, 17], "pass": [0, 1, 10, 11, 15, 17], "inject": 0, "genral": 0, "templat": 0, "all": [0, 1, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17], "ar": [0, 1, 7, 9, 10, 11, 12, 13, 15, 16, 17], "written": [0, 1, 17], "base": [0, 1, 7, 8, 9, 10, 12, 17], "directori": [0, 1, 9, 11, 13, 17], "sub": [0, 1, 15], "sh": [0, 13, 15], "an": [0, 1, 7, 9, 12, 13, 15, 16, 17], "thi": [0, 1, 7, 9, 10, 11, 12, 13, 15, 16, 17], "i": [0, 1, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17], "run": [0, 1, 8, 9, 12, 13, 15], "larg": 0, "100k": 0, "atom": [0, 1, 11, 16, 17], "supermuc": 0, "ng": 0, "bin": [0, 11, 13], "bash": [0, 11], "sbatch": 0, "job": [0, 1], "name": [0, 1, 7, 9, 11, 13, 17], "ligpair": [0, 11], "o": [0, 11], "x": [0, 1, 11, 13, 15, 16, 17], "j": [0, 11], "out": [0, 1, 11, 17], "e": [0, 1, 7, 9, 11, 16, 17], "err": 0, "d": [0, 1], "node": [0, 1, 15, 17], "130": 0, "task": 0, "per": [0, 1, 15, 17], "48": 0, "requeu": 0, "export": [0, 11, 13], "none": [0, 1, 7, 8, 10, 17], "get": [0, 1, 9, 11, 12], "env": 0, "account": 0, "xxx": 0, "partit": 0, "time": [0, 1, 9, 10, 11, 17], "10": [0, 11, 13, 17], "00": [0, 9, 13, 15, 17], "modul": [0, 13, 14], "load": [0, 13], "slurm_setup": 0, "14": [0, 1, 17], "gcc8": 0, "impi": 0, "nodes_per_namd": [0, 15], "cpus_per_namd": [0, 15], "480": 0, "echo": 0, "chang": [0, 1, 9, 11, 13, 16, 17], "line": [0, 1, 11, 12, 15], "point": [0, 1, 9, 10, 11, 17], "your": [0, 13], "project": 0, "ties_dir": [0, 11, 15], "hppf": 0, "work": [0, 1, 13, 15], "pn98ve": 0, "di67rov": 0, "test_ti": 0, "studi": [0, 17], "prot": [0, 9], "l2": 0, "l1": 0, "com": [0, 11, 13, 17], "cd": [0, 11, 13, 15, 17], "replica": [0, 1, 9, 10, 11, 15, 17], "conf": [0, 1, 15, 17], "stage": [0, 1, 15, 17], "0": [0, 1, 7, 9, 11, 13, 15, 16, 17], "do": [0, 1, 9, 10, 11, 13, 15, 16, 17], "lambda": [0, 7, 9, 10, 11, 14, 15, 16, 17], "05": [0, 15, 16, 17], "1": [0, 1, 7, 9, 11, 13, 15, 16, 17], "4": [0, 11, 13, 15, 16, 17], "5": [0, 11, 15, 16, 17], "6": [0, 11, 13, 15, 16, 17], "7": [0, 11, 13, 15, 16, 17], "8": [0, 11, 13, 15, 16, 17], "9": [0, 11, 15, 16, 17], "95": [0, 15, 16, 17], "srun": [0, 15, 17], "n": [0, 11, 13, 15, 17], "tclmain": [0, 15, 17], "sim": [0, 1, 15, 17], "sleep": [0, 15, 17], "done": [0, 11, 15, 17], "wait": [0, 11, 15, 17], "first": [0, 16, 17], "20": [0, 12, 16], "could": [0, 13, 15, 17], "adapt": 0, "smaller": [0, 15], "10k": 0, "follow": [0, 11, 13, 15, 17], "13": [0, 11, 17], "45": 0, "micro": 0, "scale": [0, 15, 16], "up": [0, 8, 10, 11, 15, 17], "simul": [0, 1, 10, 11, 12, 13, 15, 16], "bsub": [0, 11], "p": [0, 11, 13], "w": [0, 1, 11], "nnode": [0, 11], "alloc_flag": [0, 11], "gpudefault": [0, 11], "smt1": [0, 11], "test": [0, 1, 7, 13, 15, 17], "otest": 0, "etest": 0, "ls_subcwd": 0, "path": [0, 1, 8, 9, 10, 11, 13, 16], "gpf": 0, "alpin": 0, "scratch": 0, "adw62": [0, 13], "chm155": 0, "ties_test": 0, "miniconda": [0, 13], "ties_md": [0, 11, 12, 13, 15, 17], "ethan": [0, 17], "zero_sum": [0, 17], "leg1": 0, "cuda": [0, 1, 13, 15], "168": 0, "date": [0, 17], "jsrun": [0, 11, 15], "smpiarg": [0, 11, 15], "off": [0, 11, 13, 15, 17], "c": [0, 11, 13, 15], "g": [0, 1, 11, 15, 16, 17], "b": [0, 11, 13, 15, 16], "pack": [0, 11, 15], "config_fil": [0, 1, 8, 11, 15, 17], "cfg": [0, 1, 8, 11, 15, 17], "exp_nam": [0, 1, 11, 15, 17], "sys_solv": [0, 15, 17], "windows_mask": [0, 1, 11, 15, 17], "node_id": [0, 1, 11, 15, 17], "namd3": [0, 1, 15, 17], "thetagpu": [0, 15], "cobalt": 0, "A": [0, 1, 7, 11, 13, 16, 17], "t": [0, 1], "100": 0, "q": 0, "full": 0, "mpirun": 0, "lu": 0, "theta": 0, "fs0": 0, "softwar": [0, 12], "openmpi": 0, "compbioaffin": 0, "awad": 0, "namd_3": 0, "0alpha9_linux": 0, "x86_64": [0, 13], "multicor": 0, "node1": 0, "sed": 0, "1q": 0, "cobalt_nodefil": 0, "node2": 0, "2q": 0, "many_rep": 0, "mcl1": 0, "l18": 0, "l39": 0, "host": 0, "cpu": [0, 1, 13, 17], "set": [0, 1, 7, 8, 11, 13, 15, 16, 17], "bind": [0, 1, 12, 16, 17], "core": [0, 17], "np": [0, 1, 9, 10], "devic": [0, 1, 15, 17], "30": 0, "40": 0, "50": [0, 17], "60": 0, "70": 0, "80": [0, 13], "90": 0, "alchem": [0, 1, 7, 11, 12, 13, 15, 17], "window": [0, 1, 9, 10, 15, 16, 17], "us": [0, 1, 7, 9, 10, 11, 12, 13, 15, 16, 17], "onli": [0, 15, 17], "each": [0, 1, 9, 10, 11, 15, 16, 17], "addition": 0, "gpu": [0, 1, 13, 15, 17], "idl": 0, "For": [0, 15, 16, 17], "real": 0, "world": 0, "applic": [0, 16, 17], "need": [0, 1, 11, 13, 15, 16, 17], "current": [0, 1, 17], "build": [0, 1, 8, 9, 11, 12, 13, 17], "creat": [0, 1, 15, 16, 17], "python": [0, 13, 15, 17], "veri": 0, "help": [0, 1], "would": [0, 15, 17], "allow": [0, 1, 9, 13, 15, 17], "u": [0, 1, 11, 15], "import": [0, 11], "__name__": 0, "__main__": 0, "acc_nam": 0, "how": [0, 1, 15, 16, 17], "mani": [0, 1, 15, 17], "want": [0, 1, 9, 10, 11, 13, 15], "what": [0, 1, 9, 10, 11, 15, 16, 17], "thermodynam": [0, 10, 11, 12, 13, 16], "leg": [0, 1, 9, 11, 17], "mai": [0, 1, 11, 17], "have": [0, 7, 11, 13, 15, 17], "differ": [0, 1, 7, 9, 10, 11, 13, 15, 17], "wall": 0, "where": [0, 1, 8, 9, 10, 13, 17], "binari": 0, "namd3_ex": 0, "cwd": [0, 1, 11], "getcwd": [0, 11], "give": [0, 11], "lig": [0, 9, 11], "differnt": 0, "wall_tim": 0, "els": [0, 10], "open": [0, 12], "join": [0, 11], "thetagpu_": 0, "format": [0, 1, 8, 9, 17], "f": 0, "header": [0, 1, 11], "instal": [0, 12, 17], "locat": [0, 1, 9, 15, 17], "read": [0, 1, 7, 8, 9, 10, 11, 13, 17], "file": [0, 1, 8, 9, 10, 11, 15, 16, 17], "rang": [0, 1, 15], "move": [0, 17], "iter": [0, 1, 9, 11], "over": [0, 1, 9, 11, 15, 17], "minim": [0, 1], "nvt": [0, 1], "eq": 0, "npt": [0, 1], "product": [0, 1, 17], "sim0": 0, "sim1": [0, 1, 17], "sim2": 0, "sim3": 0, "count": 0, "lam": [0, 1], "siulat": 0, "rep": [0, 1, 9, 10], "number": [0, 1, 9, 10, 17], "next": [0, 11], "when": [0, 1, 11, 15, 16, 17], "fill": [0, 17], "sure": [0, 11], "between": [0, 1, 10, 11, 13, 15, 16, 17], "finish": [0, 1, 11], "eng_script": [1, 14], "cfg_script": [1, 2], "namd_sub": [1, 2], "namd_sub_split": [1, 2], "openmm_sub_split": [1, 2], "openmmtool": [1, 13, 14], "alchemi": [1, 14], "ties_analysi": [1, 14, 17], "engin": [1, 8, 17], "namd": [1, 8, 12, 16, 17], "openmm": [1, 7, 8, 11, 12, 17], "method": [1, 8, 9, 11, 12, 15, 16, 17], "fep": [1, 8, 9, 13, 17], "config": [1, 11, 14, 15, 16, 17], "class": [1, 7, 8, 9, 10, 17], "run_typ": [1, 11, 17], "period": [1, 17], "true": [1, 7, 10], "platform": [1, 13], "kwarg": [1, 7], "object": [1, 7, 8, 9, 10], "control": [1, 15, 16, 17], "protocol": [1, 11, 12, 13], "initi": [1, 11], "variabl": [1, 7], "call": [1, 11, 16, 17], "function": [1, 7, 8, 9, 10, 17], "start": [1, 12, 13, 15], "input": [1, 10, 11, 12, 16], "script": [1, 11, 12, 15, 17], "paramet": [1, 7, 8, 9, 10, 16, 17], "str": [1, 7, 8, 9], "flag": 1, "sai": 1, "should": [1, 10, 13, 15, 16, 17], "dynam": [1, 15, 17], "experi": [1, 11, 17], "complex": [1, 11, 17], "pdb": [1, 11, 17], "prmtop": [1, 17], "list": [1, 9, 10, 16, 17], "int": [1, 9, 10], "which": [1, 7, 9, 10, 11, 13, 15, 16, 17], "float": [1, 7, 9, 10, 16], "id": [1, 9, 13], "denot": [1, 17], "execut": [1, 11, 15, 17], "contain": [1, 8, 9, 10, 11, 16, 17], "end": [1, 10, 15, 16], "boolean": [1, 10], "determin": [1, 13, 16, 17], "custom": [1, 13], "schedul": [1, 9, 10, 11, 16, 17], "sting": 1, "valu": [1, 7, 9, 10, 11, 13, 15, 16, 17], "opencl": [1, 13], "dict": [1, 8, 10], "from": [1, 7, 9, 10, 11, 15, 16, 17], "properti": 1, "box_typ": [1, 17], "type": [1, 17], "box": [1, 15, 17], "being": [1, 7, 15, 16, 17], "cube": [1, 17], "truncatedoctahedron": [1, 17], "rhombicdodecahedron": [1, 17], "na": [1, 17], "manual": [1, 17], "return": [1, 7, 8, 9, 10, 13], "build_results_dir": 1, "folder": [1, 17], "helper": [1, 9], "output": [1, 9, 10, 13, 15, 17], "param": [1, 8], "elec_edg": [1, 16, 17], "electrostat": [1, 16, 17], "potenti": [1, 9, 10, 17], "begin": [1, 16, 17], "stop": [1, 17], "appear": [1, 9, 16, 17], "molecular": [1, 15, 17], "go": 1, "prefix": [1, 13, 17], "expect": [1, 17], "experiment_nam": 1, "get_opt": [1, 8, 11], "print": [1, 8, 11], "global_lambda": [1, 15, 16, 17], "global": [1, 7, 17], "reps_per_exec": [1, 11, 15, 17], "program": [1, 11, 16, 17], "parallel": [1, 12, 17], "setup": [1, 8, 12, 17], "split_run": 1, "tiesmd": [1, 11], "subset": 1, "bool": [1, 9, 10], "split": 1, "ster_edg": [1, 16, 17], "lennard_jon": [1, 17], "total_rep": [1, 11, 15, 17], "total": [1, 15, 17], "update_cfg": [1, 11], "congig": 1, "after": [1, 17], "api": [1, 11], "made": [1, 17], "write_analysis_cfg": 1, "configur": [1, 11, 13, 15, 17], "analysi": [1, 8, 9, 10, 12, 13], "write_namd_eq": 1, "eq1": [1, 17], "eq2": [1, 17], "equilibr": [1, 17], "write_namd_min": 1, "eq0": [1, 17], "write_namd_prod": 1, "write_namd_script": 1, "write_namd_submiss": 1, "submiss": [1, 11, 12, 15, 17], "hpc": [1, 11, 12, 15, 17], "archer2": 1, "write_openmm_submiss": 1, "get_box_vector": 1, "comput": [1, 10, 11, 13, 17], "vector": [1, 17], "know": 1, "edg": [1, 17], "length": [1, 17], "defin": [1, 7, 17], "bix": 1, "unit": [1, 7, 9, 10, 17], "angstrom": [1, 7, 17], "vec3": 1, "basi": 1, "get_header_and_run": 1, "namd_vers": [1, 9, 17], "prep": [1, 11, 17], "inspect": [1, 9, 11], "best": [1, 16], "guess": 1, "version": [1, 7, 9, 13, 17], "ani": [1, 9, 11, 13, 15, 17], "one": [1, 9, 10, 11, 15, 16, 17], "eah": 1, "nice_print": [1, 8], "pad": [1, 8], "alchsi": 1, "temperatur": [1, 9, 10, 17], "pressur": [1, 17], "constraint_fil": [1, 17], "constraint_column": [1, 17], "basis_vector": 1, "input_typ": [1, 17], "amber": [1, 17], "absolut": 1, "fals": [1, 7, 10], "debug": 1, "free": [1, 10, 12, 15, 16, 17], "energi": [1, 10, 12, 15, 16, 17], "calcul": [1, 9, 10, 11, 12, 13, 15, 17], "thermostat": [1, 17], "barostat": [1, 17], "detail": [1, 11, 16, 17], "constraint": [1, 11, 17], "beta": 1, "occup": [1, 17], "column": [1, 17], "mbar": [1, 10], "explicit": 1, "cell": [1, 15, 17], "pme": 1, "cutoffnonperiod": 1, "remov": [1, 9, 10], "forc": [1, 13], "nonbond": 1, "note": [1, 11, 15, 17], "gromac": 1, "add_consraint": 1, "add": [1, 17], "constrain": [1, 17], "dure": [1, 17], "construct": [1, 17], "boundari": 1, "condit": 1, "amend_original_posit": 1, "posit": [1, 17], "updat": [1, 9, 11, 13], "store": [1, 8], "clash": 1, "found": [1, 11, 16, 17], "initialis": 1, "build_simul": 1, "device_id": 1, "take": [1, 11, 15, 16, 17], "index": [1, 9, 17], "integr": [1, 10, 12, 13, 16], "debug_forc": 1, "while": 1, "maintain": 1, "ensembl": [1, 15, 17], "modifi": [1, 7, 11, 17], "get_gradi": 1, "val": 1, "context": [1, 7], "h": 1, "analitic_ster": 1, "gradient": [1, 9, 10, 17], "r": 1, "eg": 1, "lambda_electrostatics_appear": 1, "finit": [1, 16], "analyt": 1, "steric": 1, "experiment": 1, "numer": [1, 10, 15], "get_intersect_angl": 1, "appear_idx": 1, "disappear_idx": 1, "idx": 1, "angl": 1, "straddl": 1, "region": 1, "dissapear": 1, "get_intersect_bond": 1, "bond": 1, "get_intersect_tors": 1, "torsion": 1, "disappear": [1, 9, 16, 17], "rebuild_tors": 1, "intersect_tors": 1, "rebuild": 1, "without": [1, 15, 17], "non": [1, 16, 17], "physic": [1, 12, 16, 17], "fulli": [1, 16], "result": [1, 9, 10, 11, 13, 15, 17], "nan": 1, "eval": 1, "refer": [1, 13, 15, 16, 17], "set_context_to_st": 1, "param_v": 1, "see": [1, 11, 12, 15, 17], "shift_alchemical_posit": 1, "small": [1, 16, 17], "pertub": 1, "resolv": 1, "caus": [1, 11], "overlap": [1, 10], "test_sim": 1, "find": 1, "undefin": [1, 7], "pdb_line": 1, "let": 1, "address": 1, "readabl": 1, "system_id": 1, "inform": [1, 9, 10, 12, 13, 15, 16, 17], "repeat": [1, 15], "add_simulation_report": 1, "total_step": 1, "save": [1, 10, 17], "report": 1, "step": [1, 11, 13, 16], "dcd": 1, "equilibri": 1, "save_fil": 1, "perform": [1, 9, 10, 11, 13, 17], "state": [1, 7, 10, 16, 17], "indic": 1, "whether": 1, "ha": [1, 7, 13, 15, 17], "get_alchemical_atom": 1, "pdb_file": 1, "pull": 1, "temp": [1, 9, 10, 13], "factor": [1, 17], "logic": 1, "get_constraint": 1, "fuction": 1, "row": 1, "data": [1, 9, 10], "strength": 1, "relax": 1, "preproduct": [1, 17], "equili_step": 1, "equili_state_fil": 1, "meta": 1, "wrap": 1, "togeth": 1, "equilib": 1, "remove_simulation_report": 1, "strip": 1, "simulate_system": 1, "alch_si": 1, "mask": [1, 15], "niter": 1, "steps_per_it": 1, "1000": 1, "main": [1, 8, 17], "Will": 1, "assign": 1, "worker": 1, "pre": [1, 17], "collect": [1, 9, 12], "requir": [1, 7, 17], "grad": [1, 10], "disk": [1, 11], "sampl": [1, 9, 10, 12, 16, 17], "2f": 1, "much": [1, 17], "argv": 1, "entri": [1, 10, 17], "command": [1, 11, 12, 13, 15], "argument": [1, 11], "read_config": [1, 8, 11], "arg": [1, 7], "global_lamb": 1, "fraction": 1, "electr": 1, "kei": [1, 17], "lambda_": 1, "_": [1, 7], "termin": [1, 8, 11, 13], "update_attrs_from_schedul": [1, 9], "attribut": 1, "self": [1, 9, 10], "lambda_sterics_appear": [1, 9], "appear_func": 1, "evalu": 1, "two": [1, 11, 12, 15, 16, 17], "y": [1, 16, 17], "disappear_func": 1, "get_lin": 1, "global_lam": 1, "numpi": [1, 9, 10], "arrai": [1, 9, 10], "interpol": 1, "across": [1, 15], "except": 7, "alchemicalstateerror": 7, "globalparametererror": 7, "error": [7, 10, 13, 15], "rais": 7, "alchemicalst": 7, "modifiedabsolutealchemicalfactori": 7, "consistent_except": 7, "switch_width": 7, "quantiti": 7, "alchemical_pme_treat": 7, "exact": [7, 17], "alchemical_rf_treat": 7, "switch": 7, "disable_alchemical_dispersion_correct": 7, "split_alchemical_forc": 7, "absolutealchemicalfactori": 7, "super": 7, "modifiedalchemicalst": 7, "parameters_name_suffix": 7, "globalparameterst": 7, "apply_to_context": 7, "put": [7, 11], "If": [7, 11, 13, 15, 17], "doe": [7, 11, 13, 17], "apply_to_system": 7, "check_system_consist": 7, "check": 7, "It": [7, 17], "consist": 7, "classmethod": 7, "from_system": 7, "constructor": 7, "specifi": [7, 9, 10, 15, 17], "search": 7, "parameter_nam": 7, "repres": 7, "same": [7, 11, 15, 17], "get_alchemical_vari": 7, "variable_nam": 7, "warn": 7, "deprec": 7, "get_function_vari": 7, "instead": [7, 13], "variable_valu": 7, "enter": 7, "mathemat": 7, "express": 7, "alchemicalfunct": 7, "enslav": 7, "arbitrari": 7, "lambda_angl": 7, "interv": 7, "standard": [7, 10, 17], "lambda_bond": 7, "lambda_electrostat": 7, "lambda_ster": 7, "lambda_tors": 7, "set_alchemical_paramet": 7, "new_valu": 7, "given": [7, 17], "those": 7, "remain": 7, "new": [7, 9, 11], "set_alchemical_vari": 7, "set_function_vari": 7, "analysis_cfg": 8, "hold": [8, 9], "chose": 8, "win_mask": 9, "distribut": [9, 10], "rep_convg": [9, 10], "sampling_convg": [9, 10], "vdw_a": 9, "vdw_d": 9, "ele_a": 9, "ele_d": 9, "dir": 9, "writen": [9, 17], "dg": [9, 10], "individu": [9, 10, 11, 15], "intermedi": [9, 16], "you": [9, 11, 13, 17], "wish": [9, 11, 13, 17], "converg": 9, "amount": 9, "describ": 9, "vdw": 9, "elec": 9, "alch": [9, 14], "collate_data": 9, "data_root": 9, "protein": [9, 11, 12, 16, 17], "ligand": [9, 11, 12, 16, 17], "thermo": 9, "run_analysi": 9, "kelvin": [9, 10, 17], "implement": [9, 12, 15, 16], "stdev": 9, "get_it": 9, "file_loc": 9, "get_replica": 9, "sort": 9, "specif": [9, 12, 16, 17], "rep0": 9, "get_window": 9, "lambda_0": 9, "read_alch_fil": 9, "file_path": 9, "namd_ver": 9, "ver": 9, "old": 9, "look": [9, 10, 13, 15, 17], "12": [9, 11, 17], "about": [9, 10, 17], "easili": 9, "queryabl": 9, "etc": 9, "fep_combine_rep": 9, "openmm_ti": [9, 17], "combin": [9, 15], "seri": [9, 10], "concaten": 9, "mbar_analysi": 10, "analysis_dir": 10, "decorrel": 10, "averag": 10, "befor": 10, "process": 10, "u_kln": 10, "n_k": 10, "mask_window": 10, "rep_id": 10, "matrix": 10, "associ": [10, 11, 16, 17], "decorrelate_data": 10, "decorrol": 10, "turpl": 10, "plot_overlap_mat": 10, "mat": 10, "plot": [10, 16], "replica_analysi": 10, "consid": [10, 15, 16, 17], "trajectori": 10, "bootstrap": 10, "sem": [10, 11], "ti_analysi": 10, "deviat": [10, 17], "intergr": 10, "lambda_paramet": 10, "varianc": 10, "plot_du_by_dl": 10, "du": 10, "dlam": 10, "v": [10, 13], "includ": [10, 11, 16], "ci": 10, "compute_bs_error": 10, "var": 10, "boot": 10, "strap": 10, "get_lam_diff": 10, "lambda_arrai": 10, "adjac": 10, "discuss": [11, 15], "section": [11, 12, 13, 17], "outlin": [11, 12, 15, 17], "shorten": 11, "respect": [11, 17], "transform": [11, 13, 16, 17], "anoth": [11, 16], "water": 11, "now": [11, 15], "insid": 11, "\u03b4": 11, "g_": 11, "mutation1": 11, "mutation2": 11, "equal": [11, 15, 17], "\u03b4\u03b4": 11, "idea": [11, 15], "relat": 11, "shown": [11, 16], "figur": [11, 16], "more": [11, 12, 15, 16, 17], "public": [11, 16], "cournia": 11, "et": 11, "al": 11, "whole": 11, "cycl": 11, "must": [11, 13, 17], "hand": 11, "ties20": [11, 16, 17], "principl": 11, "addit": [11, 15], "howev": [11, 15], "normal": 11, "avoid": 11, "rapid": 11, "crystal": 11, "structur": [11, 17], "conform": 11, "earli": 11, "close": 11, "contact": 11, "also": [11, 13, 15, 17], "so": [11, 13, 15, 17], "expens": 11, "computation": 11, "recommend": [11, 16, 17], "design": 11, "both": [11, 17], "parameter": 11, "hybrid": [11, 16], "browser": 11, "altern": [11, 15, 17], "pair": 11, "ligand_ff_nam": 11, "gaff2": 11, "ligand_net_charg": 11, "md_engin": 11, "liganda": 11, "mol2": 11, "antechamber_prepare_mol2": 11, "ligandb": 11, "ensur": [11, 13, 15, 17], "make_atom_names_uniqu": 11, "turn": [11, 17], "superimpos": 11, "sinc": 11, "declar": 11, "prepare_input": 11, "receptor": 11, "protein_ff": 11, "leaprc": 11, "ff14sb": 11, "re": 11, "prepar": [11, 12], "built": [11, 12, 16, 17], "order": [11, 13], "via": 11, "At": 11, "map": 11, "thermodynamic_leg": [11, 17], "wa": [11, 13, 17], "seen": 11, "good": [11, 17], "default": [11, 13, 17], "quickli": 11, "cli": [11, 14], "thermo_leg": 11, "args_dict": 11, "chm155_001": 11, "120": 11, "oligpair": 11, "eligpair": 11, "expr": 11, "abov": [11, 13, 15, 17], "intern": 11, "yield": [11, 13], "3": [11, 12, 13, 15, 16, 17], "11": [11, 13], "These": [11, 13, 15, 16, 17], "onc": 11, "last": 11, "ties_ana": [11, 17], "Then": [11, 13], "our": [11, 12, 15, 17], "like": [11, 13, 15, 17], "minu": 11, "care": [11, 17], "compar": 11, "depend": [11, 17], "again": [11, 15], "dat": [11, 17], "packag": [12, 13, 14, 17], "stand": 12, "enhanc": [12, 16], "accur": 12, "reproduc": 12, "rel": [12, 16], "pleas": [12, 17], "theori": 12, "within": [12, 13, 17], "avail": [12, 13, 15, 17], "onlin": [12, 17], "sourc": [12, 17], "explain": 12, "websit": 12, "builder": 12, "linux": 12, "ppc64le": 12, "tutori": [12, 13], "bfe": 12, "background": 12, "pathwai": [12, 17], "github": [12, 13, 15, 17], "conda": 13, "manag": 13, "assum": 13, "wget": 13, "http": [13, 17], "repo": 13, "continuum": 13, "io": 13, "miniconda3": 13, "latest": 13, "chmod": 13, "match": 13, "machin": 13, "permiss": 13, "final": 13, "page": [13, 15, 17], "With": [13, 17], "forg": 13, "verifi": 13, "m": 13, "testinstal": 13, "older": 13, "simtk": 13, "instanc": [13, 15, 17], "wrong": 13, "cudatoolkit": 13, "happen": 13, "produc": [13, 17], "git": [13, 17], "revis": 13, "130124a3f9277b054ec40927360a6ad20c8f5fa6": 13, "There": [13, 17], "successfulli": 13, "cuda_error_unsupported_ptx_vers": 13, "222": 13, "median": 13, "30571e": 13, "06": 13, "76359e": 13, "05194e": 13, "07": 13, "toler": 13, "critic": [13, 15, 17], "correct": [13, 17], "replac": 13, "particular": [13, 17], "One": 13, "appropri": 13, "nvidia": 13, "smi": 13, "460": 13, "driver": 13, "persist": 13, "bu": 13, "disp": 13, "volatil": 13, "uncorr": 13, "ecc": 13, "fan": 13, "perf": 13, "pwr": 13, "usag": 13, "cap": 13, "memori": 13, "util": [13, 17], "mig": 13, "quadro": 13, "m1000m": 13, "00000000": 13, "01": 13, "50c": 13, "p5": 13, "435mib": 13, "2002mib": 13, "top": 13, "right": 13, "my": [13, 17], "correctli": [13, 17], "download": [13, 17], "clone": [13, 17], "ucl": [13, 17], "cc": [13, 17], "pymabr": 13, "therefor": [13, 17], "until": 13, "around": 13, "copi": 13, "elsewher": 13, "skip": 13, "And": 13, "mkdir": 13, "openmmtools_instal": 13, "powerpc": 13, "pip": 13, "featur": 13, "tree": 13, "subpackag": 14, "content": 14, "submodul": [14, 17], "domain": 15, "kind": 15, "spatial": 15, "were": 15, "decompos": 15, "difficult": 15, "achiev": [15, 17], "than": 15, "ones": 15, "focu": 15, "aleator": 15, "inher": 15, "chaotic": 15, "commun": 15, "other": 15, "embarrassingli": 15, "problem": [15, 17], "easi": 15, "likewis": 15, "remaind": 15, "explor": 15, "th": 15, "spit": 15, "separ": [15, 16, 17], "tell": [15, 17], "otherwis": 15, "ident": 15, "uniqu": [15, 17], "further": 15, "inclus": 15, "exclus": 15, "maximum": 15, "clariti": 15, "multipl": [15, 17], "resourc": 15, "alloc": 15, "handl": 15, "case": [15, 17], "messag": 15, "interfac": 15, "mpi": 15, "vari": [15, 17], "univers": 15, "solut": [15, 17], "independ": [15, 17], "someth": [15, 17], "notic": [15, 17], "loop": [15, 17], "65": 15, "anecdot": 15, "less": 15, "crash": 15, "135": 15, "extens": 15, "provis": 15, "hardwar": 15, "exhaust": 15, "suggest": [15, 17], "consult": [15, 16, 17], "document": [15, 16, 17], "comprehens": 15, "molecul": [16, 17], "coval": 16, "dual": 16, "topologi": [16, 17], "approach": 16, "complementari": 16, "involv": [16, 17], "chemic": 16, "moieti": 16, "destroi": 16, "stratifi": 16, "along": 16, "chosen": 16, "keep": 16, "track": 16, "\u03bb": [16, 17], "As": [16, 17], "three": [16, 17], "tune": 16, "wai": 16, "comma": [16, 17], "exactli": [16, 17], "lennard": 16, "jone": 16, "lj": 16, "interact": 16, "second": [16, 17], "graphic": 16, "guid": 16, "intend": 17, "submit": 17, "code": 17, "navig": 17, "essenti": 17, "support": 17, "them": 17, "place": 17, "understand": 17, "renam": 17, "anyth": 17, "fix": 17, "novel": 17, "desir": 17, "occasion": 17, "alongsid": 17, "valid": 17, "300": 17, "atmospher": 17, "4n": 17, "sampling_per_window": 17, "nanosecond": 17, "2n": 17, "equili_per_window": 17, "estim": 17, "least": 17, "evoc": 17, "parallelis": 17, "space": 17, "con": 17, "beta_factor": 17, "edge_length": 17, "nanomet": 17, "cell_basis_vec1": 17, "cell_basis_vec2": 17, "cell_basis_vec3": 17, "effect": 17, "alpha": 17, "simpl": 17, "imag": 17, "show": 17, "appli": 17, "alwai": 17, "leap": 17, "log": 17, "preform": 17, "popul": 17, "size": 17, "typic": 17, "possibli": 17, "drug": 17, "solvent": 17, "later": 17, "why": 17, "present": 17, "almost": 17, "readi": 17, "invok": 17, "just": 17, "taken": 17, "either": 17, "halt": 17, "advanc": 17, "below": 17, "silent": 17, "ignor": 17, "integ": 17, "request": 17, "_alpha": 17, "By": 17, "intens": 17, "pc": 17, "head": 17, "lambda_x": 17, "repi": 17, "analys": 17, "zero": 17, "sum": 17, "smallest": 17, "branch": 17, "mind": 17, "6x1": 17, "lot": 17, "calcualt": 17, "128": 17, "might": 17, "explicitli": 17, "wide": 17, "solv": 17, "issu": 17, "miss": 17, "under": 17, "exp": 17, "theoret": 17, "\u03b4g": 17, "kcal": 17, "mol": 17, "unknown": 17, "measur": 17, "left": 17, "becaus": 17, "carri": 17, "special": 17, "infer": 17, "dictionari": 17, "methodologi": 17, "openmm_fep": 17, "023": 17, "003": 17, "076": 17}, "objects": {"": [[1, 0, 0, "-", "TIES_MD"]], "TIES_MD": [[1, 0, 0, "-", "TIES"], [1, 0, 0, "-", "alch"], [1, 0, 0, "-", "cli"], [2, 0, 0, "-", "eng_scripts"], [1, 0, 0, "-", "lambdas"], [7, 0, 0, "-", "openmmtools"], [8, 0, 0, "-", "ties_analysis"]], "TIES_MD.TIES": [[1, 1, 1, "", "TIES"], [1, 4, 1, "", "get_box_vectors"], [1, 4, 1, "", "get_header_and_run"], [1, 4, 1, "", "nice_print"]], "TIES_MD.TIES.TIES": [[1, 2, 1, "", "box_type"], [1, 3, 1, "", "build_results_dirs"], [1, 2, 1, "", "elec_edges"], [1, 2, 1, "", "engine"], [1, 2, 1, "", "exp_name"], [1, 3, 1, "", "get_options"], [1, 2, 1, "", "global_lambdas"], [1, 2, 1, "", "reps_per_exec"], [1, 3, 1, "", "run"], [1, 3, 1, "", "setup"], [1, 2, 1, "", "split_run"], [1, 2, 1, "", "ster_edges"], [1, 2, 1, "", "total_reps"], [1, 3, 1, "", "update_cfg"], [1, 3, 1, "", "write_analysis_cfg"], [1, 3, 1, "", "write_namd_eq"], [1, 3, 1, "", "write_namd_min"], [1, 3, 1, "", "write_namd_prod"], [1, 3, 1, "", "write_namd_scripts"], [1, 3, 1, "", "write_namd_submissions"], [1, 3, 1, "", "write_openmm_submission"]], "TIES_MD.alch": [[1, 1, 1, "", "AlchSys"], [1, 1, 1, "", "PDB_line"], [1, 1, 1, "", "System_ID"], [1, 4, 1, "", "add_simulation_reporters"], [1, 4, 1, "", "equilibriation"], [1, 4, 1, "", "get_alchemical_atoms"], [1, 4, 1, "", "get_constraints"], [1, 4, 1, "", "minimization"], [1, 4, 1, "", "preproduction"], [1, 4, 1, "", "remove_simulation_reporters"], [1, 4, 1, "", "simulate_system"]], "TIES_MD.alch.AlchSys": [[1, 3, 1, "", "add_consraints"], [1, 3, 1, "", "amend_original_positions"], [1, 3, 1, "", "build_simulation"], [1, 3, 1, "", "debug_force"], [1, 3, 1, "", "get_gradients"], [1, 3, 1, "", "get_intersect_angles"], [1, 3, 1, "", "get_intersect_bonds"], [1, 3, 1, "", "get_intersect_torsions"], [1, 3, 1, "", "rebuild_torsion"], [1, 3, 1, "", "set_context_to_state"], [1, 3, 1, "", "shift_alchemical_positions"], [1, 3, 1, "", "test_sim"]], "TIES_MD.cli": [[1, 4, 1, "", "main"], [1, 4, 1, "", "read_config"]], "TIES_MD.eng_scripts": [[3, 0, 0, "-", "cfg_scripts"], [4, 0, 0, "-", "namd_sub"], [5, 0, 0, "-", "namd_sub_split"], [6, 0, 0, "-", "openmm_sub_split"]], "TIES_MD.lambdas": [[1, 1, 1, "", "Lambdas"], [1, 4, 1, "", "appear_func"], [1, 4, 1, "", "disappear_func"], [1, 4, 1, "", "get_line"]], "TIES_MD.lambdas.Lambdas": [[1, 3, 1, "", "update_attrs_from_schedule"]], "TIES_MD.openmmtools": [[7, 0, 0, "-", "alchemy"]], "TIES_MD.openmmtools.alchemy": [[7, 5, 1, "", "AlchemicalStateError"], [7, 1, 1, "", "ModifiedAbsoluteAlchemicalFactory"], [7, 1, 1, "", "ModifiedAlchemicalState"]], "TIES_MD.openmmtools.alchemy.ModifiedAlchemicalState": [[7, 3, 1, "", "apply_to_context"], [7, 3, 1, "", "apply_to_system"], [7, 3, 1, "", "check_system_consistency"], [7, 3, 1, "", "from_system"], [7, 3, 1, "", "get_alchemical_variable"], [7, 3, 1, "", "get_function_variable"], [7, 6, 1, "", "lambda_angles"], [7, 6, 1, "", "lambda_bonds"], [7, 6, 1, "", "lambda_electrostatics"], [7, 6, 1, "", "lambda_sterics"], [7, 6, 1, "", "lambda_torsions"], [7, 3, 1, "", "set_alchemical_parameters"], [7, 3, 1, "", "set_alchemical_variable"], [7, 3, 1, "", "set_function_variable"]], "TIES_MD.ties_analysis": [[8, 0, 0, "-", "config"], [9, 0, 0, "-", "engines"], [10, 0, 0, "-", "methods"], [8, 0, 0, "-", "ties_analysis"]], "TIES_MD.ties_analysis.config": [[8, 1, 1, "", "Config"], [8, 4, 1, "", "read_config"]], "TIES_MD.ties_analysis.config.Config": [[8, 3, 1, "", "get_options"]], "TIES_MD.ties_analysis.engines": [[9, 0, 0, "-", "namd"], [9, 0, 0, "-", "openmm"]], "TIES_MD.ties_analysis.engines.namd": [[9, 1, 1, "", "NAMD"], [9, 4, 1, "", "get_iter"], [9, 4, 1, "", "get_replica"], [9, 4, 1, "", "get_window"], [9, 4, 1, "", "read_alch_file"]], "TIES_MD.ties_analysis.engines.namd.NAMD": [[9, 3, 1, "", "collate_data"], [9, 3, 1, "", "run_analysis"]], "TIES_MD.ties_analysis.engines.openmm": [[9, 1, 1, "", "Lambdas"], [9, 1, 1, "", "OpenMM"], [9, 4, 1, "", "get_replica"], [9, 4, 1, "", "get_window"]], "TIES_MD.ties_analysis.engines.openmm.Lambdas": [[9, 3, 1, "", "update_attrs_from_schedule"]], "TIES_MD.ties_analysis.engines.openmm.OpenMM": [[9, 3, 1, "", "collate_data"], [9, 3, 1, "", "run_analysis"]], "TIES_MD.ties_analysis.methods": [[10, 0, 0, "-", "FEP"], [10, 0, 0, "-", "TI"]], "TIES_MD.ties_analysis.methods.FEP": [[10, 1, 1, "", "MBAR_Analysis"]], "TIES_MD.ties_analysis.methods.FEP.MBAR_Analysis": [[10, 3, 1, "", "analysis"], [10, 3, 1, "", "decorrelate_data"], [10, 3, 1, "", "plot_overlap_mat"], [10, 3, 1, "", "replica_analysis"]], "TIES_MD.ties_analysis.methods.TI": [[10, 1, 1, "", "TI_Analysis"], [10, 4, 1, "", "compute_bs_error"], [10, 4, 1, "", "get_lam_diff"]], "TIES_MD.ties_analysis.methods.TI.TI_Analysis": [[10, 3, 1, "", "analysis"], [10, 3, 1, "", "intergrate"], [10, 3, 1, "", "plot_du_by_dl"]], "TIES_MD.ties_analysis.ties_analysis": [[8, 1, 1, "", "Analysis"], [8, 4, 1, "", "main"], [8, 4, 1, "", "nice_print"]], "TIES_MD.ties_analysis.ties_analysis.Analysis": [[8, 3, 1, "", "run"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:property", "3": "py:method", "4": "py:function", "5": "py:exception", "6": "py:attribute"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "property", "Python property"], "3": ["py", "method", "Python method"], "4": ["py", "function", "Python function"], "5": ["py", "exception", "Python exception"], "6": ["py", "attribute", "Python attribute"]}, "titleterms": {"hpc": 0, "submiss": 0, "script": 0, "namd": [0, 9, 15], "openmm": [0, 9, 13, 15], "3": 0, "ties_md": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14], "packag": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "subpackag": [1, 2, 8], "submodul": [1, 7, 8, 9, 10], "ti": [1, 10, 12, 13, 15], "modul": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "alch": 1, "cli": 1, "lambda": 1, "content": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12], "eng_script": [2, 3, 4, 5, 6], "cfg_script": 3, "namd_sub": 4, "namd_sub_split": 5, "openmm_sub_split": 6, "openmmtool": 7, "alchemi": 7, "ties_analysi": [8, 9, 10], "config": 8, "engin": 9, "method": 10, "fep": 10, "bind": 11, "free": 11, "energi": 11, "tutori": [11, 17], "gener": 11, "bfe": 11, "background": 11, "setup": 11, "run": [11, 17], "analysi": [11, 17], "welcom": 12, "md": [12, 13], "document": 12, "code": 12, "instal": 13, "linux": 13, "ppc64le": 13, "parallel": 15, "theori": 16, "outlin": 16, "alchem": 16, "calcul": 16, "pathwai": 16, "get": 17, "start": 17, "input": 17, "command": 17, "line": 17, "simul": 17, "prepar": 17}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})