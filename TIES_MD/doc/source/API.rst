TIES MD API
===========

``TIES MD`` can be used on the command line but for greater automation we also provide an API that exposes some options that
may be routinely changed during setup.

API
---

Here we detail all the options in the API and what should be passed. The options that were previously on the command line
can be passed into the TIES class like so::

    from TIES_MD import TIES
    import os
    md = TIES(cwd=os.path.abspath('./my_ties_sims'), windows_mask=[0,1], rep_id=0, exp_name='sys_solv')

Once the TIES class is constructed the options that were previously in TIES.cfg can now be set as attributes of the TIES
class like so::

    # openmm.unit is needed to set values with units
    from openmm import unit

    #string for the molecular dynamics engine (openmm/namd2.14/namd3)
    md.engine = 'openmm'

    #Target temperature for the thermostat
    md.temperature = 300.0*unit.kelvin

    #Target pressure for barostat
    md.pressure = 1.0*unit.atmosphere

    #How much production sampling to run per alchemical window (4ns recommended)
    md.sampling_per_window = 0.04*unit.nanosecond

    #How much equilibration to run per alchemical window (2ns recommended)
    md.equili_per_window = 0.002*unit.nanosecond

    #List for which estimators to use.
    md.methods = ['FEP', 'TI']

    #How many total replicas of each window are run (we recommend at least 5).
    md.total_reps = 3

    #Boolean for if we will split all replicas into separate runs. (True for maximum parallelism)
    md.split_run = False

    #List for where in lambda schedule (0->1) should the electrostatic potentials begin, stop appearing.
    md.elec_edges = [0.5, 1.0]

    #List for where in lambda schedule (0->1) should the Lennard_Jones potentials begin, stop appearing.
    md.ster_edges = [0.0, 0.5]

    #List for the value the global controlling parameter takes in each window.
    md.global_lambdas = [0.00, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.00]

    #String for name of the pdb file with constraints in the build directory, i.e. 'cons.pdb'
    md.constraint_file = None

    #String for column in pdb are constraints provided valid options are 'occupancy'/'beta_factor'.
    md.constraint_column = None

    #String for what input type is provided, only AMBER supported.
    md.input_type = 'AMBER'

    #list of x, y, z floats for box vectors of this simulation, unit Angstrom.
    md.cell_basis_vec1 = [34.55, 0.0, 0.0]
    md.cell_basis_vec2 = [-11.516722937414105, 32.574214501232206, 0.0]
    md.cell_basis_vec3 = [-11.516722937414105, -16.287105279373797, 28.21009840448772]


Finally there are three additional options that don't appear in TIES.cfg which are ``sub_header``,
``pre_run_line`` and ``run_line``. These three option can be used as follows::

    #A header to a the submission script that will be written for this job.
    md.sub_header = """#Example script for Summit OpenMM
    #BSUB -P CHM155_001
    #BSUB -W 120
    #BSUB -nnodes 13
    #BSUB -alloc_flags "gpudefault smt1"
    #BSUB -J LIGPAIR
    #BSUB -o oLIGPAIR.%J
    #BSUB -e eLIGPAIR.%J"""

    #The prefix for the run line of this job.
    md.pre_run_line = 'jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 '

    #the ties_md run line of this job.
    md.run_line = 'ties_md --config_file=$ties_dir/TIES.cfg --windows_mask=$lambda,$(expr $lambda + 1) --node_id=$i'

The use of these three setting will produce a submission script which looks like::

    #!/bin/bash
    #Example script for Summit OpenMM
    #BSUB -P CHM155_001
    #BSUB -W 120
    #BSUB -nnodes 13
    #BSUB -alloc_flags "gpudefault smt1"
    #BSUB -J LIGPAIR
    #BSUB -o oLIGPAIR.%J
    #BSUB -e eLIGPAIR.%J

    export ties_dir="ties/ties-ligandA-ligandB/lig"
    cd $ties_dir

    for lambda in 0 1 2 3 4 5 6 7 8 9 10 11 12; do
      for i in {0..5}; do
            jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --windows_mask=$lambda,$(expr $lambda + 1) --node_id=$i&
        done
        done
    wait

If ``sub_header``, ``pre_run_line`` and ``run_line`` are not set ``TIES_MD`` will make a best guess for a submission script.
Ideally only small modification should be need to run using the best guess scripts. Any tweaks that are applied to get the
scripts working can then be passed into ``sub_header``, ``pre_run_line`` and ``run_line`` for future system setups. For
general ideas on how to make submission scripts see :ref:`HPC Submission scripts`.

