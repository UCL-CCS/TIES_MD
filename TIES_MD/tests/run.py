#!/usr/bin/env python

import pytest, tempfile, os, site, shutil

def main():
    with tempfile.TemporaryDirectory() as tmpdirname:
        print('created temporary directory', tmpdirname)
        test_dir = os.path.join(site.getsitepackages()[0], 'TIES_MD', 'tests')
        tmp_test_dir = os.path.join(tmpdirname, 'TIES_MD', 'unit_testing', 'tests')
        shutil.copytree(test_dir, tmp_test_dir)
        pytest.main([os.path.join(tmp_test_dir)])

if __name__ == '__main__':
    main()



