quick:
	g++ -shared -std=c++2a -undefined dynamic_lookup `python3 -m pybind11 --includes` fold.cpp -o fold`python3-config --extension-suffix` -D PYTHON
