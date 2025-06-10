#!/bin/bash

# Test script for MS²Rescore with Percolator Docker image

echo "🧪 Testing MS²Rescore with Percolator Docker image..."
echo "=================================================="

# Docker image tag
IMAGE_TAG="registry.gitlab.com/bzhanglab/ms2rescore:3.1.5"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to run test and report result
run_test() {
    local test_name="$1"
    local command="$2"
    echo -e "${YELLOW}Testing: ${test_name}${NC}"
    
    if eval "$command"; then
        echo -e "${GREEN}✅ PASS: ${test_name}${NC}"
        echo ""
        return 0
    else
        echo -e "${RED}❌ FAIL: ${test_name}${NC}"
        echo ""
        return 1
    fi
}

# Test 1: Check if Docker image exists
run_test "Docker image exists" \
    "docker images | grep -q 'registry.gitlab.com/bzhanglab/ms2rescore'"

# Test 2: Basic container startup
run_test "Container starts successfully" \
    "docker run --rm ${IMAGE_TAG} --help > /dev/null 2>&1"

# Test 3: MS²Rescore version
echo -e "${YELLOW}Testing: MS²Rescore version${NC}"
docker run --rm ${IMAGE_TAG} --version
echo ""

# Test 4: Percolator installation
run_test "Percolator is accessible" \
    "docker run --rm --entrypoint='' ${IMAGE_TAG} bash -c 'percolator --help > /dev/null 2>&1'"

# Test 5: Percolator version
echo -e "${YELLOW}Testing: Percolator version${NC}"
docker run --rm --entrypoint="" ${IMAGE_TAG} bash -c "percolator --help | head -5"
echo ""

# Test 6: Python dependencies
run_test "Python dependencies installed" \
    "docker run --rm --entrypoint='' ${IMAGE_TAG} python -c 'import ms2rescore, ms2pip, deeplc, mokapot; print(\"All imports successful\")'"

# Test 7: Feature generators available
echo -e "${YELLOW}Testing: Available feature generators${NC}"
docker run --rm --entrypoint="" ${IMAGE_TAG} python -c "
import ms2rescore.feature_generators as fg
generators = [name for name in dir(fg) if not name.startswith('_')]
print('Available generators:', generators)
"
echo ""

# Test 8: Rescoring engines available
echo -e "${YELLOW}Testing: Available rescoring engines${NC}"
docker run --rm --entrypoint="" ${IMAGE_TAG} python -c "
import ms2rescore.rescoring_engines as re
engines = [name for name in dir(re) if not name.startswith('_')]
print('Available engines:', engines)
"
echo ""

# Test 9: Test configuration parsing
echo -e "${YELLOW}Testing: Configuration file parsing${NC}"
docker run --rm --entrypoint="" -v $(pwd):/test ${IMAGE_TAG} python -c "
import sys
sys.path.append('/test')
try:
    from ms2rescore.core import read_config
    config = read_config('/test/percolator-config.toml')
    print('✅ Configuration parsed successfully')
    print('Rescoring engine:', list(config.get('rescoring_engine', {}).keys()))
except Exception as e:
    print('❌ Configuration parsing failed:', str(e))
    sys.exit(1)
"
echo ""

# Test 10: Full help output
echo -e "${YELLOW}Testing: Full help output${NC}"
docker run --rm ${IMAGE_TAG} --help
echo ""

echo "🎉 Testing completed!"
echo "If all tests passed, your Docker image is ready to use!" 