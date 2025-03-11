all: build test

build:
	cmake -S. -Bbuild
	cmake --build build

test: build
	cd build && \
	ctest -V -j 10