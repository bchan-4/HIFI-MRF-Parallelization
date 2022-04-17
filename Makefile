HIFIserial: src/HIFI_serial.cpp src/HIFI_options.cpp
	g++ -O4 -o src/HIFIserial src/HIFI_serial.cpp src/HIFI_options.cpp -lm
HIFI: src/HIFI.cpp src/HIFI_options.cpp
	mpicc -O3 -o src/HIFI src/HIFI.cpp src/HIFI_options.cpp -lm -fopenmp
HIFIAdvanced: src/HIFI_Advanced.cpp src/HIFI_options.cpp
	mpicc -O3 -o src/HIFIAdvanced src/HIFI_Advanced.cpp src/HIFI_options.cpp -lm -fopenmp
callPeaks: src/callPeaks.cpp
	g++ -O4 -o src/callPeaks src/callPeaks.cpp -lm