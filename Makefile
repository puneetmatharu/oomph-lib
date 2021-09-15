NPROCS:=1
OS := $(shell uname -s)

ifeq ($(OS),Linux)
	NPROCS := $(shell nproc)
endif
ifeq ($(OS),Darwin)
	NPROCS := $(shell sysctl -n hw.logicalcpu)
endif


setup:
	echo "Can't use this yet. WIP. Need to make sure it runs safely."
	mkdir build && cd build
	if command -v ninja &> /dev/null; then
		cmake -G Ninja ..
		ninja && ninja install
	else
		cmake ..
		make && make install
	fi

check:
	touch validation.log && rm -f validation.log
	cd demo_drivers
	mkdir build
	cd build
	if command -v ninja &> /dev/null; then
		cmake -G Ninja ..
	else
		cmake ..
	fi
	ctest -j $(NPROCS) --timeout 10000 --no-label-summary

cmake-format:
	set -f
	find . \( -name '*.cmake' -o -name 'CMakeLists.txt' \) -exec cmake-format -c .cmake-format.py -i {} \;
	set +f

clang-format:
	echo "Now starting clang-format run on src/ and demo_drivers/ folder..."
	./scripts/run_clang_format_on_all.sh ./src
	echo "Finished running clang-format."