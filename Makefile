regular:
	g++ -shared -std=c++2a -undefined_dynamic_lookup `python3 -m pybind11 --includes` fold.cpp -o fold`python3-config --extension-suffix` -D PYTHON -fPIC
quick:
	g++ -shared -std=c++2a -O3 -undefined_dynamic_lookup `python3 -m pybind11 --includes` fold.cpp -o fold`python3-config --extension-suffix` -D PYTHON -fPIC 
