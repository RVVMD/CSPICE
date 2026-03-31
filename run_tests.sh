#!/bin/bash
# CSPICE Test Runner Script
# Builds and runs all tests, reports results

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
TESTS_DIR="${SCRIPT_DIR}/tests"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}===========================================${NC}"
echo -e "${BLUE}   CSPICE Test Runner${NC}"
echo -e "${BLUE}===========================================${NC}"
echo ""

# Build the project with tests
echo -e "${YELLOW}Building CSPICE with tests...${NC}"
cd "${BUILD_DIR}"
cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON
make -j4

echo ""
echo -e "${YELLOW}Running unit tests...${NC}"
echo ""

# Run the test executable
if [ -f "${BUILD_DIR}/bin/cspice_tests" ]; then
    "${BUILD_DIR}/bin/cspice_tests"
    TEST_EXIT_CODE=$?
else
    echo -e "${RED}Error: Test executable not found${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}Running netlist simulations...${NC}"
echo ""

# Run transformer test
if [ -f "${SCRIPT_DIR}/transformer_test.cir" ]; then
    echo -e "${BLUE}Running transformer_test.cir...${NC}"
    "${BUILD_DIR}/bin/cspice" "${SCRIPT_DIR}/transformer_test.cir"
    echo ""
fi

# Run three-phase test
if [ -f "${SCRIPT_DIR}/three_phase_inrush.cir" ]; then
    echo -e "${BLUE}Running three_phase_inrush.cir...${NC}"
    "${BUILD_DIR}/bin/cspice" "${SCRIPT_DIR}/three_phase_inrush.cir"
    echo ""
fi

echo -e "${BLUE}===========================================${NC}"
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
else
    echo -e "${RED}Some tests failed!${NC}"
fi
echo -e "${BLUE}===========================================${NC}"

exit $TEST_EXIT_CODE
