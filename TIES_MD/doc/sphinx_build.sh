pip install ../..
sphinx-apidoc -f -o ./source/ ..
sphinx-apidoc .. -o ./source/ -f 
make html
