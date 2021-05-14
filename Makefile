regular:
	g++ -shared -std=c++17 -undefined_dynamic_lookup `python3 -m pybind11 --includes` fold.cpp -o fold`python3-config --extension-suffix` -D PYTHON -fPIC
quick:
	g++ -shared -std=c++17 -O3 -undefined_dynamic_lookup `python3 -m pybind11 --includes` fold.cpp -o fold`python3-config --extension-suffix` -D PYTHON -fPIC
cpp:
	g++ -std=c++17 fold.cpp -g -o fold
