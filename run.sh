#!/usr/bin/env bash
# Build (if needed) and launch the Botany realtime viewer.
#
# First run configures CMake and downloads a handful of C++ libraries via
# FetchContent (glm, GLFW, glad, imgui, Catch2, evolve). Takes a few minutes.
# Subsequent runs are incremental and usually under a second to launch.
#
# Any arguments are passed through to botany_realtime, e.g.
#     ./run.sh --color auxin
#     ./run.sh --load-plant saves/plant_42_tick1000.tree

set -euo pipefail

cd "$(dirname "$0")"

if ! command -v cmake >/dev/null 2>&1; then
    echo "error: cmake not found." >&2
    echo "  macOS:  brew install cmake" >&2
    echo "  Debian: sudo apt install cmake build-essential" >&2
    exit 1
fi

if [ ! -f build/CMakeCache.txt ]; then
    cmake -B build
fi
cmake --build build -j

./build/botany_realtime "$@"
