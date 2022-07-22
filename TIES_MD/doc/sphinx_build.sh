pip install ../..
sphinx-apidoc –implicit-namespaces -f -o ./source/ ..
sphinx-apidoc .. -o ./source/ -f 
make html
