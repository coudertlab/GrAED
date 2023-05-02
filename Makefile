all: src/ads_grid.cpp
	c++ -std=c++11 -I./include -O3 src/ads_grid.cpp -o graed
clean:
	rm -f graed
