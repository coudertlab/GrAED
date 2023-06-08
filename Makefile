all: src/ads_grid.cpp
	c++ -std=c++11 -I./include -O3 src/ads_grid.cpp -o graed
	c++ -std=c++11 -I./include -O3 src/grid.cpp -o pgrid
clean:
	rm -f graed
	rm -f pgrid
