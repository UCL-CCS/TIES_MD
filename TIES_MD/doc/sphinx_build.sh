pip install ../..
sphinx-apidoc -f -o ./source/ .. ../eng_scripts/namd_many_rep/*, ../eng_scripts/openmm/*, ../eng_scripts/namd_single_rep/*, ../eng_scripts/*, ../eng_scripts, ../openmmtools/*, ../openmmtools
sphinx-apidoc .. -o ./source/ -f
make html
