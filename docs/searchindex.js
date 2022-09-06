Search.setIndex({"docnames": ["API", "HPC_submissions", "TIES_MD", "TIES_MD.eng_scripts", "TIES_MD.eng_scripts.cfg_scripts", "TIES_MD.eng_scripts.namd_sub", "TIES_MD.eng_scripts.namd_sub_split", "TIES_MD.eng_scripts.openmm_sub", "TIES_MD.eng_scripts.openmm_sub_split", "TIES_MD.openmmtools", "TIES_MD.ties_analysis", "TIES_MD.ties_analysis.engines", "TIES_MD.ties_analysis.methods", "binding_free_energies", "index", "install", "modules", "parallelization", "theory", "tutorial"], "filenames": ["API.rst", "HPC_submissions.rst", "TIES_MD.rst", "TIES_MD.eng_scripts.rst", "TIES_MD.eng_scripts.cfg_scripts.rst", "TIES_MD.eng_scripts.namd_sub.rst", "TIES_MD.eng_scripts.namd_sub_split.rst", "TIES_MD.eng_scripts.openmm_sub.rst", "TIES_MD.eng_scripts.openmm_sub_split.rst", "TIES_MD.openmmtools.rst", "TIES_MD.ties_analysis.rst", "TIES_MD.ties_analysis.engines.rst", "TIES_MD.ties_analysis.methods.rst", "binding_free_energies.rst", "index.rst", "install.rst", "modules.rst", "parallelization.rst", "theory.rst", "tutorial.rst"], "titles": ["TIES MD API", "HPC Submission scripts", "TIES_MD package", "TIES_MD.eng_scripts package", "TIES_MD.eng_scripts.cfg_scripts package", "TIES_MD.eng_scripts.namd_sub package", "TIES_MD.eng_scripts.namd_sub_split package", "TIES_MD.eng_scripts.openmm_sub package", "TIES_MD.eng_scripts.openmm_sub_split package", "TIES_MD.openmmtools package", "TIES_MD.ties_analysis package", "TIES_MD.ties_analysis.engines package", "TIES_MD.ties_analysis.methods package", "Binding Free Energy Tutorial", "Welcome to the TIES MD documentation.", "Installation", "TIES_MD", "Parallelization", "Theory", "Tutorial"], "terms": {"can": [0, 1, 2, 9, 12, 13, 14, 15, 17, 18, 19], "us": [0, 1, 2, 9, 11, 12, 13, 14, 15, 17, 18, 19], "command": [0, 2, 13, 14, 15, 17], "line": [0, 1, 2, 13, 14, 17], "greater": 0, "autom": 0, "we": [0, 1, 2, 10, 11, 12, 13, 17, 18, 19], "also": [0, 13, 15, 17, 19], "provid": [0, 1, 2, 15, 17, 18, 19], "an": [0, 1, 2, 9, 11, 14, 15, 17, 18, 19], "expos": [0, 1], "some": [0, 1, 2, 13, 15, 17, 18, 19], "option": [0, 1, 2, 9, 10, 11, 13, 17, 18, 19], "mai": [0, 1, 2, 13, 19], "routin": 0, "chang": [0, 1, 2, 11, 13, 15, 18, 19], "dure": [0, 2, 19], "setup": [0, 2, 10, 14, 19], "here": [0, 1, 2, 13, 15, 17, 18, 19], "detail": [0, 2, 13, 18, 19], "all": [0, 1, 2, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19], "what": [0, 1, 2, 11, 12, 13, 17, 18, 19], "should": [0, 2, 12, 15, 17, 18, 19], "pass": [0, 1, 2, 12, 13, 17, 19], "The": [0, 1, 2, 9, 13, 15, 17, 18, 19], "were": [0, 17], "previous": 0, "class": [0, 2, 9, 10, 11, 12, 19], "A": [0, 1, 2, 9, 13, 15, 18, 19], "minim": [0, 1, 2, 19], "evoc": 0, "would": [0, 1, 17, 19], "look": [0, 11, 12, 15, 17, 19], "like": [0, 13, 15, 17, 19], "from": [0, 2, 9, 11, 12, 13, 17, 18, 19], "ties_md": [0, 1, 13, 14, 15, 17, 19], "import": [0, 1, 13], "o": [0, 1, 13], "cwd": [0, 1, 2], "path": [0, 1, 2, 10, 11, 12, 13, 15, 18, 19], "abspath": 0, "my_ties_sim": 0, "point": [0, 1, 2, 11, 12, 13, 19], "where": [0, 1, 2, 10, 11, 12, 15, 19], "cfg": [0, 1, 2, 10, 13, 17, 19], "file": [0, 1, 2, 10, 11, 12, 13, 17, 18, 19], "i": [0, 1, 2, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19], "locat": [0, 1, 2, 11, 19], "simul": [0, 1, 2, 12, 13, 14, 15, 17, 18], "ar": [0, 1, 2, 9, 11, 12, 13, 14, 15, 17, 18, 19], "perform": [0, 2, 11, 12, 13, 15, 19], "note": [0, 2, 13, 17, 19], "thi": [0, 1, 2, 9, 11, 12, 13, 14, 15, 17, 18, 19], "must": [0, 13, 15, 19], "absolut": [0, 2], "other": [0, 17], "argument": [0, 2, 13], "set": [0, 1, 2, 9, 10, 13, 15, 17, 18, 19], "default": [0, 13, 15, 19], "valu": [0, 2, 9, 11, 12, 13, 15, 18, 19], "If": [0, 9, 13, 15, 17, 19], "want": [0, 1, 2, 11, 12, 13, 15, 17], "code": [0, 19], "might": [0, 19], "exp_nam": [0, 1, 2, 17, 19], "sys_solv": [0, 1, 17, 19], "windows_mask": [0, 1, 2, 13, 17, 19], "0": [0, 1, 2, 9, 11, 13, 15, 17, 18, 19], "devic": [0, 1, 2, 17, 19], "rep_id": [0, 1, 2, 12, 17, 19], "These": [0, 13, 15, 17, 18, 19], "have": [0, 1, 9, 13, 15, 17, 19], "same": [0, 9, 13, 17, 19], "definit": 0, "subsect": 0, "tutori": [0, 14, 15], "onc": [0, 13], "construct": [0, 2, 19], "now": [0, 13, 17], "attribut": [0, 2], "so": [0, 13, 15, 17, 19], "openmm": [0, 2, 9, 10, 13, 14, 19], "unit": [0, 2, 9, 11, 12, 19], "need": [0, 1, 2, 13, 15, 17, 18, 19], "string": [0, 1, 2, 10, 11, 12, 19], "molecular": [0, 2, 17, 19], "dynam": [0, 2, 17, 19], "engin": [0, 2, 10, 19], "namd2": [0, 1, 2, 17, 19], "14": [0, 1, 2, 19], "namd3": [0, 1, 2, 17, 19], "target": [0, 1, 2, 15, 18, 19], "temperatur": [0, 2, 11, 12, 19], "thermostat": [0, 2, 19], "300": [0, 19], "kelvin": [0, 11, 12, 19], "pressur": [0, 2, 19], "barostat": [0, 2, 19], "1": [0, 1, 2, 9, 11, 13, 15, 17, 18, 19], "atmospher": [0, 19], "how": [0, 1, 2, 13, 17, 18, 19], "much": [0, 2, 19], "product": [0, 1, 2, 19], "sampl": [0, 2, 11, 12, 14, 18, 19], "run": [0, 1, 2, 10, 11, 14, 15, 17], "per": [0, 1, 2, 17, 19], "alchem": [0, 1, 2, 9, 13, 14, 15, 17, 19], "window": [0, 1, 2, 11, 12, 17, 18, 19], "4n": [0, 19], "recommend": [0, 13, 18, 19], "sampling_per_window": [0, 19], "04": 0, "nanosecond": [0, 19], "equilibr": [0, 2, 19], "2n": [0, 19], "equili_per_window": [0, 19], "002": 0, "list": [0, 2, 11, 12, 18, 19], "which": [0, 2, 9, 11, 12, 13, 15, 17, 18, 19], "estim": [0, 19], "method": [0, 2, 10, 11, 13, 14, 17, 18, 19], "fep": [0, 2, 10, 11, 15, 19], "mani": [0, 1, 2, 17, 19], "total": [0, 2, 17, 19], "replica": [0, 1, 2, 11, 12, 13, 17, 19], "each": [0, 1, 2, 11, 12, 13, 17, 18, 19], "least": [0, 19], "5": [0, 1, 17, 18, 19], "total_rep": [0, 13, 17, 19], "3": [0, 14, 15, 17, 18, 19], "boolean": [0, 2, 12, 19], "split": [0, 2, 19], "separ": [0, 17, 18, 19], "true": [0, 2, 9, 10, 12, 13], "maximum": [0, 17, 19], "parallel": [0, 2, 14, 19], "split_run": [0, 2, 13, 17, 19], "fals": [0, 2, 9, 12], "lambda": [0, 1, 9, 11, 12, 13, 16, 17, 18, 19], "schedul": [0, 2, 11, 12, 13, 18, 19], "electrostat": [0, 2, 18, 19], "potenti": [0, 2, 11, 12, 19], "begin": [0, 2, 18, 19], "stop": [0, 2, 19], "appear": [0, 2, 11, 18, 19], "elec_edg": [0, 2, 18, 19], "lennard_jon": [0, 2, 19], "ster_edg": [0, 2, 18, 19], "global": [0, 2, 9, 19], "control": [0, 2, 17, 18, 19], "paramet": [0, 2, 9, 10, 11, 12, 18, 19], "take": [0, 2, 13, 17, 18, 19], "global_lambda": [0, 2, 17, 18, 19], "00": [0, 1, 11, 15, 17, 19], "05": [0, 1, 17, 18, 19], "2": [0, 1, 11, 15, 17, 18, 19], "4": [0, 1, 15, 17, 18, 19], "6": [0, 1, 13, 15, 17, 18, 19], "7": [0, 1, 15, 17, 18, 19], "8": [0, 1, 15, 17, 18, 19], "9": [0, 1, 17, 18, 19], "95": [0, 1, 17, 18, 19], "name": [0, 1, 2, 9, 11, 13, 15, 19], "pdb": [0, 2, 13, 19], "constraint": [0, 2, 13, 19], "build": [0, 1, 2, 10, 11, 13, 14, 15, 19], "directori": [0, 1, 2, 11, 13, 15, 19], "e": [0, 1, 2, 9, 11, 13, 18, 19], "con": [0, 19], "constraint_fil": [0, 2, 19], "none": [0, 1, 2, 9, 10, 12, 19], "column": [0, 2, 19], "valid": [0, 19], "occup": [0, 2, 19], "beta_factor": [0, 19], "constraint_column": [0, 2, 19], "input": [0, 2, 12, 13, 14, 18], "type": [0, 2, 19], "onli": [0, 1, 17, 19], "amber": [0, 2, 19], "support": [0, 19], "input_typ": [0, 2, 19], "x": [0, 1, 2, 13, 15, 17, 18, 19], "y": [0, 2, 18, 19], "z": [0, 2], "float": [0, 2, 9, 11, 12, 18, 19], "box": [0, 2, 17, 19], "vector": [0, 2, 19], "angstrom": [0, 2, 9, 19], "cell_basis_vec1": [0, 2, 19], "34": 0, "55": 0, "cell_basis_vec2": [0, 2, 19], "11": [0, 15], "516722937414105": 0, "32": 0, "574214501232206": 0, "cell_basis_vec3": [0, 2, 19], "16": 0, "287105279373797": 0, "28": 0, "21009840448772": 0, "final": [0, 15], "three": [0, 18, 19], "addit": [0, 13, 17], "don": 0, "t": [0, 1, 2], "sub_head": [0, 1, 13], "pre_run_lin": [0, 1, 13], "run_lin": [0, 1, 13], "follow": [0, 1, 13, 15, 17, 19], "header": [0, 1, 2, 13], "submiss": [0, 2, 13, 14, 17, 19], "script": [0, 2, 13, 14, 17, 19], "written": [0, 1, 2, 19], "job": [0, 1, 2], "exampl": [0, 1, 2, 13, 15, 17, 18, 19], "summit": [0, 1, 2, 13, 17], "bsub": [0, 1, 13], "p": [0, 1, 13, 15], "chm155_001": [0, 13], "w": [0, 1, 2, 13], "120": [0, 13], "nnode": [0, 1, 13], "13": [0, 1, 13, 19], "alloc_flag": [0, 1, 13], "gpudefault": [0, 1, 13], "smt1": [0, 1, 13], "j": [0, 1, 13], "ligpair": [0, 1, 13], "oligpair": [0, 13], "eligpair": [0, 13], "prefix": [0, 15, 19], "jsrun": [0, 1, 13, 17], "smpiarg": [0, 1, 13, 17], "off": [0, 1, 13, 15, 17, 19], "n": [0, 1, 13, 15, 17, 19], "c": [0, 1, 13, 15, 17], "g": [0, 1, 2, 13, 17, 18, 19], "b": [0, 1, 13, 15, 17, 18], "pack": [0, 1, 13, 17], "config_fil": [0, 1, 2, 10, 13, 17, 19], "ties_dir": [0, 1, 13, 17], "expr": [0, 13], "node_id": [0, 13], "produc": [0, 15, 19], "bin": [0, 1, 15], "bash": [0, 1], "export": [0, 1, 15], "liganda": [0, 13], "ligandb": [0, 13], "lig": [0, 1, 11, 13], "cd": [0, 1, 13, 15, 19], "10": [0, 1, 15], "12": [0, 11], "do": [0, 1, 2, 11, 12, 13, 15, 17, 18, 19], "done": [0, 1, 17, 19], "wait": [0, 1, 17, 19], "make": [0, 1, 2, 12, 13, 17], "best": [0, 2, 18], "guess": [0, 2], "ideal": 0, "small": [0, 2, 18, 19], "modif": 0, "ani": [0, 2, 11, 13, 15, 17, 19], "tweak": 0, "appli": [0, 19], "get": [0, 1, 2, 11, 13, 14], "work": [0, 1, 2, 15, 17], "futur": 0, "system": [0, 1, 2, 9, 13, 14, 15, 17, 18, 19], "For": [0, 1, 17, 18, 19], "gener": [0, 1, 10, 11, 14, 17, 18, 19], "idea": [0, 13, 17], "see": [0, 2, 13, 14, 17, 19], "hpc": [0, 2, 13, 14, 17, 19], "variou": 1, "ti": [1, 10, 11, 13, 16, 18, 19], "md": [1, 13, 17, 18, 19], "attempt": [1, 15], "automat": [1, 13, 19], "write": [1, 2, 13, 15, 17], "sensibl": 1, "archer": [1, 17], "In": [1, 13, 15, 17, 19], "user": [1, 2, 13, 15, 18, 19], "own": [1, 13, 19], "whichev": 1, "cluster": [1, 13], "thei": [1, 2, 13], "prefer": 1, "To": [1, 13, 15, 17, 18, 19], "aid": 1, "api": [1, 2, 13, 14], "call": [1, 2, 13, 18, 19], "inject": 1, "templat": 1, "base": [1, 2, 9, 10, 11, 12, 14, 19], "sub": [1, 2, 17], "sh": [1, 15, 17], "larg": 1, "100k": 1, "atom": [1, 2, 13, 18, 19], "supermuc": 1, "ng": 1, "sbatch": 1, "out": [1, 2, 13, 19], "err": 1, "d": [1, 2], "node": [1, 2, 17, 19], "130": 1, "task": 1, "48": 1, "requeu": 1, "env": 1, "account": 1, "xxx": 1, "partit": 1, "time": [1, 2, 11, 12, 13, 19], "modul": [1, 15, 16], "load": [1, 15], "slurm_setup": 1, "gcc8": 1, "impi": 1, "nodes_per_namd": [1, 17], "cpus_per_namd": [1, 17], "480": 1, "echo": 1, "your": [1, 15], "project": 1, "hppf": 1, "pn98ve": 1, "di67rov": 1, "test_ti": 1, "studi": [1, 19], "prot": [1, 11], "l2": 1, "l1": 1, "com": [1, 13, 15, 19], "conf": [1, 2, 17, 19], "stage": [1, 2, 17, 19], "srun": [1, 17, 19], "tclmain": [1, 17, 19], "sleep": [1, 17, 19], "first": [1, 18, 19], "20": [1, 14, 18, 19], "could": [1, 15, 17, 19], "adapt": 1, "smaller": [1, 17], "10k": 1, "45": 1, "micro": 1, "scale": [1, 17, 18], "up": [1, 10, 12, 13, 17, 19], "test": [1, 2, 9, 15, 17, 19], "otest": 1, "etest": 1, "ls_subcwd": 1, "gpf": 1, "alpin": 1, "scratch": 1, "adw62": [1, 15], "chm155": 1, "ties_test": 1, "miniconda": [1, 15], "ethan": [1, 19], "zero_sum": [1, 19], "leg1": 1, "cuda": [1, 2, 15, 17], "168": 1, "date": [1, 19], "thetagpu": [1, 17], "cobalt": 1, "100": 1, "q": 1, "full": 1, "mpirun": 1, "lu": 1, "theta": 1, "fs0": 1, "softwar": [1, 14], "openmpi": 1, "compbioaffin": 1, "awad": 1, "namd_3": 1, "0alpha9_linux": 1, "x86_64": [1, 15], "multicor": 1, "node1": 1, "sed": 1, "1q": 1, "cobalt_nodefil": 1, "node2": 1, "2q": 1, "many_rep": 1, "mcl1": 1, "l18": 1, "l39": 1, "host": 1, "cpu": [1, 2, 15, 19], "bind": [1, 2, 14, 18, 19], "core": [1, 19], "np": [1, 2, 11, 12], "30": 1, "40": 1, "50": [1, 19], "60": 1, "70": 1, "80": [1, 15], "90": 1, "addition": 1, "gpu": [1, 2, 15, 17, 19], "idl": 1, "real": 1, "world": 1, "applic": [1, 18, 19], "current": [1, 2, 19], "creat": [1, 2, 17, 18, 19], "python": [1, 15, 17, 19], "veri": 1, "help": [1, 2], "allow": [1, 2, 11, 15, 17, 19], "u": [1, 2, 13, 17], "__name__": 1, "__main__": 1, "acc_nam": 1, "thermodynam": [1, 12, 13, 14, 15, 18], "leg": [1, 2, 11, 13, 19], "differ": [1, 2, 9, 11, 12, 13, 15, 17, 19], "wall": 1, "binari": 1, "namd3_ex": 1, "getcwd": [1, 13], "give": [1, 13], "wall_tim": 1, "els": [1, 12], "open": [1, 14], "join": [1, 13], "thetagpu_": 1, "format": [1, 2, 10, 11, 13, 19], "f": 1, "instal": [1, 13, 14, 19], "read": [1, 2, 9, 10, 11, 12, 13, 15, 19], "rang": [1, 2, 17], "move": [1, 19], "iter": [1, 2, 11, 13], "over": [1, 2, 11, 13, 17, 19], "nvt": [1, 2, 19], "eq": 1, "npt": [1, 2, 19], "run0": 1, "run1": 1, "run2": 1, "run3": 1, "count": 1, "lam": [1, 2], "rep": [1, 2, 11, 12], "2f": [1, 2], "number": [1, 2, 11, 12, 19], "next": [1, 13], "when": [1, 2, 13, 17, 18, 19], "fill": [1, 19], "sure": 1, "between": [1, 2, 12, 13, 15, 17, 18, 19], "sim": [1, 2], "finish": [1, 2, 13, 19], "eng_script": [2, 16], "cfg_script": [2, 3], "namd_sub": [2, 3], "namd_sub_split": [2, 3], "openmm_sub": [2, 3], "openmm_sub_split": [2, 3], "openmmtool": [2, 15, 16], "alchemi": [2, 16], "ties_analysi": [2, 16, 19], "namd": [2, 10, 14, 18, 19], "config": [2, 13, 16, 17, 18, 19], "complex": [2, 13, 19], "run_typ": [2, 13, 19], "period": 2, "platform": [2, 15], "kwarg": [2, 9], "object": [2, 9, 10, 11, 12], "protocol": [2, 13, 14, 15], "initi": [2, 13], "variabl": [2, 9], "function": [2, 9, 10, 11, 12, 19], "start": [2, 14, 15, 17], "str": [2, 9, 10, 11], "experi": [2, 13, 19], "prmtop": [2, 19], "flag": 2, "sai": 2, "int": [2, 11, 12, 19], "id": [2, 11, 15], "denot": [2, 19], "execut": [2, 13, 17, 19], "contain": [2, 10, 11, 12, 13, 18, 19], "end": [2, 12, 17, 18], "determin": [2, 15, 17, 18, 19], "custom": [2, 15], "sting": 2, "opencl": [2, 15], "dict": [2, 10, 12], "build_results_dir": 2, "folder": [2, 19], "helper": [2, 11], "output": [2, 11, 12, 13, 15, 19], "param": [2, 10], "return": [2, 9, 10, 11, 12, 15], "properti": 2, "1st": 2, "basi": 2, "cell": [2, 17], "compon": 2, "2nd": 2, "3rd": 2, "go": 2, "get_opt": [2, 10, 13], "print": [2, 10, 13], "tiesmd": [2, 13], "subset": 2, "bool": [2, 11, 12], "update_cfg": 2, "after": [2, 19], "made": [2, 19], "write_analysis_cfg": 2, "configur": [2, 13, 15, 17, 19], "analysi": [2, 10, 11, 12, 14, 15], "write_namd_eq": 2, "eq1": 2, "eq2": 2, "write_namd_min": 2, "eq0": 2, "write_namd_prod": 2, "sim1": 2, "write_namd_script": 2, "write_namd_submiss": 2, "archer2": 2, "write_openmm_submiss": 2, "get_box_vector": 2, "box_typ": 2, "comput": [2, 12, 13, 15, 19], "know": 2, "edg": 2, "length": 2, "defin": [2, 9, 19], "bix": 2, "vec3": 2, "get_header_and_run": 2, "namd_vers": [2, 11], "num_window": 2, "prep": [2, 13, 19], "inspect": [2, 11, 13], "version": [2, 9, 11, 15], "one": [2, 11, 12, 13, 17, 18, 19], "nice_print": [2, 10], "pad": [2, 10], "read_config": [2, 10], "disk": 2, "arg": [2, 9], "alchsi": 2, "basis_vector": 2, "debug": 2, "free": [2, 12, 14, 17, 18, 19], "energi": [2, 12, 14, 17, 18, 19], "calcul": [2, 11, 12, 13, 14, 15, 17, 19], "beta": 2, "mbar": [2, 12], "explicit": 2, "pme": 2, "cutoffnonperiod": 2, "remov": [2, 11, 12], "forc": [2, 15], "nonbond": 2, "gromac": 2, "add_consraint": 2, "add": [2, 19], "constrain": [2, 19], "boundari": 2, "condit": 2, "amend_original_posit": 2, "posit": [2, 19], "updat": [2, 11, 15], "store": [2, 10], "clash": 2, "found": [2, 13, 18, 19], "initialis": 2, "build_simul": 2, "device_id": 2, "index": [2, 11, 19], "integr": [2, 12, 14, 15, 18], "debug_forc": 2, "while": 2, "maintain": 2, "ensembl": [2, 17, 19], "modifi": [2, 9, 13, 19], "get_gradi": 2, "val": 2, "context": [2, 9], "h": 2, "analitic_ster": 2, "gradient": [2, 11, 12, 19], "r": 2, "eg": 2, "lambda_electrostatics_appear": 2, "finit": [2, 18], "analyt": 2, "steric": 2, "experiment": 2, "numer": [2, 12, 17], "get_intersect_angl": 2, "appear_idx": 2, "disappear_idx": 2, "idx": 2, "angl": 2, "straddl": 2, "region": 2, "dissapear": 2, "get_intersect_bond": 2, "bond": 2, "get_intersect_tors": 2, "torsion": 2, "disappear": [2, 11, 18, 19], "rebuild_tors": 2, "intersect_tors": 2, "rebuild": 2, "without": [2, 17, 19], "non": [2, 18, 19], "physic": [2, 14, 18, 19], "fulli": [2, 18], "result": [2, 11, 12, 13, 15, 17, 19], "nan": 2, "eval": 2, "refer": [2, 15, 17, 18, 19], "set_context_to_st": 2, "param_v": 2, "shift_alchemical_posit": 2, "pertub": 2, "resolv": 2, "caus": [2, 13], "overlap": [2, 12], "test_sim": 2, "find": 2, "undefin": [2, 9], "pdb_line": 2, "let": 2, "address": 2, "readabl": 2, "system_id": 2, "inform": [2, 11, 12, 14, 15, 17, 18, 19], "repeat": [2, 17], "add_simulation_report": 2, "total_step": 2, "save": [2, 12, 19], "report": 2, "step": [2, 13, 15, 18], "dcd": 2, "equilibri": 2, "save_fil": 2, "state": [2, 9, 12, 18, 19], "indic": 2, "whether": 2, "ha": [2, 9, 15, 17, 19], "get_alchemical_atom": 2, "pdb_file": 2, "pull": 2, "temp": [2, 11, 12, 15], "factor": [2, 19], "logic": 2, "get_constraint": 2, "fuction": 2, "row": 2, "data": [2, 11, 12], "strength": 2, "relax": 2, "preproduct": [2, 19], "equili_step": 2, "equili_state_fil": 2, "meta": 2, "wrap": 2, "togeth": 2, "equilib": 2, "remove_simulation_report": 2, "strip": 2, "simulate_system": 2, "alch_si": 2, "mask": [2, 17], "niter": 2, "steps_per_it": 2, "1000": 2, "main": [2, 10, 19], "Will": 2, "assign": 2, "worker": 2, "pre": [2, 19], "collect": [2, 11, 14], "requir": [2, 9, 19], "grad": [2, 12], "argv": 2, "entri": [2, 12, 19], "program": [2, 13, 18, 19], "global_lamb": 2, "fraction": 2, "electr": 2, "kei": [2, 19], "lambda_": 2, "_": [2, 9], "termin": [2, 10, 13, 15], "update_attrs_from_schedul": [2, 11], "self": [2, 11, 12], "lambda_sterics_appear": [2, 11], "appear_func": 2, "evalu": 2, "two": [2, 13, 14, 17, 18, 19], "disappear_func": 2, "get_lin": 2, "global_lam": 2, "numpi": [2, 11, 12], "arrai": [2, 11, 12], "interpol": 2, "across": [2, 17], "except": 9, "alchemicalstateerror": 9, "globalparametererror": 9, "error": [9, 12, 13, 15, 17], "rais": 9, "alchemicalst": 9, "modifiedabsolutealchemicalfactori": 9, "consistent_except": 9, "switch_width": 9, "quantiti": 9, "alchemical_pme_treat": 9, "exact": [9, 19], "alchemical_rf_treat": 9, "switch": 9, "disable_alchemical_dispersion_correct": 9, "split_alchemical_forc": 9, "absolutealchemicalfactori": 9, "super": 9, "modifiedalchemicalst": 9, "parameters_name_suffix": 9, "globalparameterst": 9, "apply_to_context": 9, "put": 9, "doe": [9, 13, 15, 19], "apply_to_system": 9, "check_system_consist": 9, "check": 9, "It": [9, 19], "consist": 9, "classmethod": 9, "from_system": 9, "constructor": 9, "specifi": [9, 11, 12, 17, 19], "search": 9, "parameter_nam": 9, "repres": 9, "get_alchemical_vari": 9, "variable_nam": 9, "warn": 9, "deprec": 9, "get_function_vari": 9, "instead": [9, 15], "variable_valu": 9, "enter": 9, "mathemat": 9, "express": 9, "alchemicalfunct": 9, "enslav": 9, "arbitrari": 9, "lambda_angl": 9, "interv": 9, "standard": [9, 12, 13, 19], "lambda_bond": 9, "lambda_electrostat": 9, "lambda_ster": 9, "lambda_tors": 9, "set_alchemical_paramet": 9, "new_valu": 9, "given": [9, 19], "those": 9, "being": [9, 17, 18, 19], "remain": 9, "new": [9, 11, 13], "set_alchemical_vari": 9, "set_function_vari": 9, "analysis_cfg": 10, "hold": [10, 11], "chose": 10, "make_exp": 10, "verbos": 10, "win_mask": 11, "distribut": [11, 12], "rep_convg": [11, 12], "sampling_convg": [11, 12], "vdw_a": 11, "vdw_d": 11, "ele_a": 11, "ele_d": 11, "dir": 11, "writen": [11, 19], "dg": [11, 12], "individu": [11, 12, 13, 17], "intermedi": [11, 18], "you": [11, 13, 15, 19], "wish": [11, 13, 15, 19], "converg": 11, "amount": 11, "describ": 11, "vdw": 11, "elec": 11, "alch": [11, 16], "collate_data": 11, "data_root": 11, "protein": [11, 13, 14, 18, 19], "ligand": [11, 13, 14, 18, 19], "thermo": 11, "run_analysi": 11, "implement": [11, 14, 17, 18], "stdev": 11, "get_it": 11, "file_loc": 11, "get_replica": 11, "sort": 11, "specif": [11, 13, 14, 18, 19], "rep0": 11, "get_window": 11, "lambda_0": 11, "read_alch_fil": 11, "file_path": 11, "namd_ver": 11, "ver": 11, "old": 11, "about": [11, 12, 19], "easili": 11, "queryabl": 11, "etc": 11, "fep_combine_rep": 11, "openmm_ti": [11, 19], "combin": [11, 17], "seri": [11, 12], "concaten": 11, "mbar_analysi": 12, "analysis_dir": 12, "decorrel": 12, "averag": 12, "befor": 12, "process": 12, "u_kln": 12, "n_k": 12, "mask_window": 12, "matrix": 12, "associ": [12, 13, 18, 19], "decorrelate_data": 12, "decorrol": 12, "turpl": 12, "plot_overlap_mat": 12, "mat": 12, "plot": [12, 18], "replica_analysi": 12, "consid": [12, 17, 18, 19], "trajectori": 12, "bootstrap": 12, "sem": [12, 13], "ti_analysi": 12, "deviat": [12, 19], "intergr": 12, "lambda_paramet": 12, "varianc": 12, "plot_du_by_dl": 12, "du": 12, "dlam": 12, "v": [12, 15], "includ": [12, 13, 18], "ci": 12, "compute_bs_error": 12, "var": 12, "boot": 12, "strap": 12, "get_lam_diff": 12, "lambda_arrai": 12, "adjac": 12, "discuss": [13, 17], "section": [13, 14, 15, 19], "outlin": [13, 14, 17, 19], "shorten": 13, "respect": [13, 19], "transform": [13, 15, 18, 19], "anoth": [13, 18], "water": 13, "insid": 13, "\u03b4": 13, "g_": 13, "mutation1": 13, "mutation2": 13, "equal": [13, 17, 19], "\u03b4\u03b4": 13, "relat": 13, "shown": [13, 18, 19], "figur": [13, 18], "more": [13, 14, 17, 18, 19], "public": [13, 18], "cournia": 13, "et": 13, "al": 13, "whole": 13, "cycl": 13, "hand": 13, "ties20": [13, 18, 19], "principl": 13, "howev": [13, 17], "normal": 13, "avoid": 13, "rapid": 13, "crystal": 13, "structur": [13, 19], "conform": 13, "earli": 13, "close": 13, "contact": 13, "expens": 13, "computation": 13, "design": 13, "both": [13, 19], "parameter": 13, "hybrid": [13, 18], "browser": 13, "altern": [13, 17, 19], "order": [13, 15], "local": [13, 19], "pleas": [13, 14, 19], "git": [13, 15, 19], "With": [13, 15, 19], "pair": 13, "ligand_ff_nam": 13, "gaff2": 13, "ligand_net_charg": 13, "md_engin": 13, "mol2": 13, "antechamber_prepare_mol2": 13, "ensur": [13, 15, 19], "make_atom_names_uniqu": 13, "turn": [13, 19], "superimpos": 13, "sinc": 13, "declar": 13, "prepare_input": 13, "receptor": 13, "protein_ff": 13, "leaprc": 13, "ff14sb": 13, "re": 13, "prepar": [13, 14], "built": [13, 14, 18, 19], "via": [13, 19], "At": 13, "map": 13, "thermodynamic_leg": [13, 19], "wa": [13, 15, 19], "seen": 13, "good": [13, 19], "quickli": 13, "thermo_leg": 13, "element": 13, "part": 13, "abov": [13, 15, 17, 19], "intern": 13, "last": 13, "ties_ana": [13, 19], "Then": [13, 15], "our": [13, 14, 17, 19], "minu": 13, "care": 13, "compar": 13, "depend": [13, 19], "again": 13, "dat": [13, 19], "mean": 13, "packag": [14, 15, 16, 19], "stand": 14, "enhanc": [14, 18], "accur": 14, "reproduc": 14, "rel": [14, 18], "theori": 14, "within": [14, 15, 19], "avail": [14, 15, 17, 19], "onlin": [14, 19], "sourc": [14, 19], "explain": 14, "websit": 14, "builder": 14, "linux": 14, "ppc64le": 14, "bfe": 14, "background": 14, "pathwai": [14, 19], "github": [14, 15, 17, 19], "conda": 15, "manag": 15, "assum": 15, "wget": 15, "http": [15, 19], "repo": 15, "continuum": 15, "io": 15, "miniconda3": 15, "latest": 15, "chmod": 15, "match": 15, "machin": 15, "permiss": 15, "page": [15, 17, 19], "forg": 15, "verifi": 15, "m": 15, "testinstal": 15, "older": 15, "simtk": 15, "instanc": [15, 17, 19], "wrong": 15, "cudatoolkit": 15, "happen": 15, "revis": 15, "130124a3f9277b054ec40927360a6ad20c8f5fa6": 15, "There": [15, 19], "successfulli": 15, "cuda_error_unsupported_ptx_vers": 15, "222": 15, "median": 15, "30571e": 15, "06": 15, "76359e": 15, "05194e": 15, "07": 15, "toler": 15, "critic": [15, 17], "correct": [15, 19], "replac": 15, "particular": [15, 19], "One": 15, "appropri": 15, "nvidia": 15, "smi": 15, "yield": 15, "460": 15, "driver": 15, "persist": 15, "bu": 15, "disp": 15, "volatil": 15, "uncorr": 15, "ecc": 15, "fan": 15, "perf": 15, "pwr": 15, "usag": 15, "cap": 15, "memori": 15, "util": 15, "mig": 15, "quadro": 15, "m1000m": 15, "00000000": 15, "01": 15, "50c": 15, "p5": 15, "435mib": 15, "2002mib": 15, "top": 15, "right": 15, "my": [15, 19], "correctli": 15, "download": [15, 19], "clone": [15, 19], "ucl": [15, 19], "cc": [15, 19], "pymabr": 15, "therefor": [15, 19], "until": 15, "around": 15, "copi": 15, "elsewher": 15, "skip": 15, "And": 15, "mkdir": 15, "openmmtools_instal": 15, "powerpc": 15, "pip": 15, "featur": 15, "tree": 15, "subpackag": 16, "content": 16, "submodul": [16, 19], "cli": 16, "domain": 17, "kind": 17, "spatial": 17, "decompos": 17, "difficult": 17, "achiev": [17, 19], "than": 17, "ones": 17, "focu": 17, "aleator": 17, "inher": 17, "chaotic": 17, "commun": 17, "embarrassingli": 17, "problem": [17, 19], "easi": 17, "likewis": 17, "remaind": 17, "explor": 17, "th": 17, "spit": 17, "tell": [17, 19], "further": 17, "inclus": 17, "exclus": 17, "sing": 17, "clariti": 17, "multipl": [17, 19], "resourc": 17, "alloc": 17, "handl": 17, "case": [17, 19], "messag": 17, "interfac": 17, "mpi": 17, "vari": [17, 19], "univers": 17, "solut": [17, 19], "independ": [17, 19], "someth": [17, 19], "notic": [17, 19], "loop": [17, 19], "65": 17, "anecdot": 17, "less": 17, "crash": 17, "135": 17, "extens": 17, "provis": 17, "hardwar": 17, "exhaust": 17, "suggest": [17, 19], "consult": [17, 18, 19], "document": [17, 18, 19], "comprehens": 17, "molecul": [18, 19], "coval": 18, "dual": 18, "topologi": [18, 19], "approach": 18, "complementari": 18, "involv": [18, 19], "chemic": 18, "moieti": 18, "destroi": 18, "stratifi": 18, "along": 18, "chosen": 18, "keep": 18, "track": 18, "\u03bb": [18, 19], "As": [18, 19], "tune": 18, "wai": 18, "comma": [18, 19], "exactli": [18, 19], "lennard": 18, "jone": 18, "lj": 18, "interact": 18, "second": [18, 19], "graphic": 18, "guid": 18, "expect": 19, "navig": 19, "essenti": 19, "place": 19, "understand": 19, "renam": 19, "anyth": 19, "fix": 19, "novel": 19, "desir": 19, "servic": 19, "later": 19, "occasion": 19, "alongsid": 19, "space": 19, "na": 19, "manual": 19, "simpl": 19, "imag": 19, "show": 19, "alwai": 19, "popul": 19, "size": 19, "typic": 19, "possibli": 19, "drug": 19, "solvent": 19, "why": 19, "present": 19, "almost": 19, "readi": 19, "invok": 19, "just": 19, "taken": 19, "either": 19, "halt": 19, "advanc": 19, "integ": 19, "By": 19, "below": 19, "silent": 19, "ignor": 19, "request": 19, "uniqu": 19, "intens": 19, "pc": 19, "head": 19, "lambda_x": 19, "repi": 19, "analys": 19, "zero": 19, "sum": 19, "smallest": 19, "branch": 19, "mind": 19, "6x1": 19, "lot": 19, "128": 19, "explicitli": 19, "wide": 19, "solv": 19, "issu": 19, "submit": 19, "miss": 19, "under": 19, "exp": 19, "theoret": 19, "\u03b4g": 19, "kcal": 19, "mol": 19, "unknown": 19, "measur": 19, "left": 19, "becaus": 19, "carri": 19, "special": 19, "infer": 19, "dictionari": 19, "methodologi": 19, "openmm_fep": 19, "023": 19, "003": 19, "076": 19}, "objects": {"": [[2, 0, 0, "-", "TIES_MD"]], "TIES_MD": [[2, 0, 0, "-", "TIES"], [2, 0, 0, "-", "alch"], [2, 0, 0, "-", "cli"], [3, 0, 0, "-", "eng_scripts"], [2, 0, 0, "-", "lambdas"], [9, 0, 0, "-", "openmmtools"], [10, 0, 0, "-", "ties_analysis"]], "TIES_MD.TIES": [[2, 1, 1, "", "TIES"], [2, 4, 1, "", "get_box_vectors"], [2, 4, 1, "", "get_header_and_run"], [2, 4, 1, "", "nice_print"], [2, 4, 1, "", "read_config"]], "TIES_MD.TIES.TIES": [[2, 2, 1, "", "build_results_dirs"], [2, 3, 1, "", "cell_basis_vec1"], [2, 3, 1, "", "cell_basis_vec2"], [2, 3, 1, "", "cell_basis_vec3"], [2, 3, 1, "", "elec_edges"], [2, 3, 1, "", "engine"], [2, 2, 1, "", "get_options"], [2, 3, 1, "", "global_lambdas"], [2, 2, 1, "", "run"], [2, 2, 1, "", "setup"], [2, 3, 1, "", "split_run"], [2, 3, 1, "", "ster_edges"], [2, 2, 1, "", "update_cfg"], [2, 2, 1, "", "write_analysis_cfg"], [2, 2, 1, "", "write_namd_eq"], [2, 2, 1, "", "write_namd_min"], [2, 2, 1, "", "write_namd_prod"], [2, 2, 1, "", "write_namd_scripts"], [2, 2, 1, "", "write_namd_submissions"], [2, 2, 1, "", "write_openmm_submission"]], "TIES_MD.alch": [[2, 1, 1, "", "AlchSys"], [2, 1, 1, "", "PDB_line"], [2, 1, 1, "", "System_ID"], [2, 4, 1, "", "add_simulation_reporters"], [2, 4, 1, "", "equilibriation"], [2, 4, 1, "", "get_alchemical_atoms"], [2, 4, 1, "", "get_constraints"], [2, 4, 1, "", "minimization"], [2, 4, 1, "", "preproduction"], [2, 4, 1, "", "remove_simulation_reporters"], [2, 4, 1, "", "simulate_system"]], "TIES_MD.alch.AlchSys": [[2, 2, 1, "", "add_consraints"], [2, 2, 1, "", "amend_original_positions"], [2, 2, 1, "", "build_simulation"], [2, 2, 1, "", "debug_force"], [2, 2, 1, "", "get_gradients"], [2, 2, 1, "", "get_intersect_angles"], [2, 2, 1, "", "get_intersect_bonds"], [2, 2, 1, "", "get_intersect_torsions"], [2, 2, 1, "", "rebuild_torsion"], [2, 2, 1, "", "set_context_to_state"], [2, 2, 1, "", "shift_alchemical_positions"], [2, 2, 1, "", "test_sim"]], "TIES_MD.cli": [[2, 4, 1, "", "main"]], "TIES_MD.eng_scripts": [[4, 0, 0, "-", "cfg_scripts"], [5, 0, 0, "-", "namd_sub"], [6, 0, 0, "-", "namd_sub_split"], [7, 0, 0, "-", "openmm_sub"], [8, 0, 0, "-", "openmm_sub_split"]], "TIES_MD.lambdas": [[2, 1, 1, "", "Lambdas"], [2, 4, 1, "", "appear_func"], [2, 4, 1, "", "disappear_func"], [2, 4, 1, "", "get_line"]], "TIES_MD.lambdas.Lambdas": [[2, 2, 1, "", "update_attrs_from_schedule"]], "TIES_MD.openmmtools": [[9, 0, 0, "-", "alchemy"]], "TIES_MD.openmmtools.alchemy": [[9, 5, 1, "", "AlchemicalStateError"], [9, 1, 1, "", "ModifiedAbsoluteAlchemicalFactory"], [9, 1, 1, "", "ModifiedAlchemicalState"]], "TIES_MD.openmmtools.alchemy.ModifiedAlchemicalState": [[9, 2, 1, "", "apply_to_context"], [9, 2, 1, "", "apply_to_system"], [9, 2, 1, "", "check_system_consistency"], [9, 2, 1, "", "from_system"], [9, 2, 1, "", "get_alchemical_variable"], [9, 2, 1, "", "get_function_variable"], [9, 6, 1, "", "lambda_angles"], [9, 6, 1, "", "lambda_bonds"], [9, 6, 1, "", "lambda_electrostatics"], [9, 6, 1, "", "lambda_sterics"], [9, 6, 1, "", "lambda_torsions"], [9, 2, 1, "", "set_alchemical_parameters"], [9, 2, 1, "", "set_alchemical_variable"], [9, 2, 1, "", "set_function_variable"]], "TIES_MD.ties_analysis": [[10, 0, 0, "-", "config"], [11, 0, 0, "-", "engines"], [12, 0, 0, "-", "methods"], [10, 0, 0, "-", "ties_analysis"]], "TIES_MD.ties_analysis.config": [[10, 1, 1, "", "Config"], [10, 4, 1, "", "read_config"]], "TIES_MD.ties_analysis.config.Config": [[10, 2, 1, "", "get_options"]], "TIES_MD.ties_analysis.engines": [[11, 0, 0, "-", "namd"], [11, 0, 0, "-", "openmm"]], "TIES_MD.ties_analysis.engines.namd": [[11, 1, 1, "", "NAMD"], [11, 4, 1, "", "get_iter"], [11, 4, 1, "", "get_replica"], [11, 4, 1, "", "get_window"], [11, 4, 1, "", "read_alch_file"]], "TIES_MD.ties_analysis.engines.namd.NAMD": [[11, 2, 1, "", "collate_data"], [11, 2, 1, "", "run_analysis"]], "TIES_MD.ties_analysis.engines.openmm": [[11, 1, 1, "", "Lambdas"], [11, 1, 1, "", "OpenMM"], [11, 4, 1, "", "get_replica"], [11, 4, 1, "", "get_window"]], "TIES_MD.ties_analysis.engines.openmm.Lambdas": [[11, 2, 1, "", "update_attrs_from_schedule"]], "TIES_MD.ties_analysis.engines.openmm.OpenMM": [[11, 2, 1, "", "collate_data"], [11, 2, 1, "", "run_analysis"]], "TIES_MD.ties_analysis.methods": [[12, 0, 0, "-", "FEP"], [12, 0, 0, "-", "TI"]], "TIES_MD.ties_analysis.methods.FEP": [[12, 1, 1, "", "MBAR_Analysis"]], "TIES_MD.ties_analysis.methods.FEP.MBAR_Analysis": [[12, 2, 1, "", "analysis"], [12, 2, 1, "", "decorrelate_data"], [12, 2, 1, "", "plot_overlap_mat"], [12, 2, 1, "", "replica_analysis"]], "TIES_MD.ties_analysis.methods.TI": [[12, 1, 1, "", "TI_Analysis"], [12, 4, 1, "", "compute_bs_error"], [12, 4, 1, "", "get_lam_diff"]], "TIES_MD.ties_analysis.methods.TI.TI_Analysis": [[12, 2, 1, "", "analysis"], [12, 2, 1, "", "intergrate"], [12, 2, 1, "", "plot_du_by_dl"]], "TIES_MD.ties_analysis.ties_analysis": [[10, 1, 1, "", "Analysis"], [10, 4, 1, "", "main"], [10, 4, 1, "", "make_exp"], [10, 4, 1, "", "nice_print"]], "TIES_MD.ties_analysis.ties_analysis.Analysis": [[10, 2, 1, "", "run"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method", "3": "py:property", "4": "py:function", "5": "py:exception", "6": "py:attribute"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"], "3": ["py", "property", "Python property"], "4": ["py", "function", "Python function"], "5": ["py", "exception", "Python exception"], "6": ["py", "attribute", "Python attribute"]}, "titleterms": {"ti": [0, 2, 12, 14, 15, 17], "md": [0, 14, 15], "api": 0, "hpc": 1, "submiss": 1, "script": 1, "namd": [1, 11, 17], "openmm": [1, 11, 15, 17], "3": 1, "ties_md": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16], "packag": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "subpackag": [2, 3, 10], "submodul": [2, 9, 10, 11, 12], "modul": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "alch": 2, "cli": 2, "lambda": 2, "content": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14], "eng_script": [3, 4, 5, 6, 7, 8], "cfg_script": 4, "namd_sub": 5, "namd_sub_split": 6, "openmm_sub": 7, "openmm_sub_split": 8, "openmmtool": 9, "alchemi": 9, "ties_analysi": [10, 11, 12], "config": 10, "engin": 11, "method": 12, "fep": 12, "bind": 13, "free": 13, "energi": 13, "tutori": [13, 19], "gener": 13, "bfe": 13, "background": 13, "setup": 13, "run": [13, 19], "analysi": [13, 19], "welcom": 14, "document": 14, "code": 14, "instal": 15, "linux": 15, "ppc64le": 15, "parallel": 17, "theori": 18, "outlin": 18, "alchem": 18, "calcul": 18, "pathwai": 18, "get": 19, "start": 19, "input": 19, "command": 19, "line": 19, "simul": 19, "prepar": 19}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})