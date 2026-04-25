#!/bin/bash
set -e

echo "------------------------------------------------"
echo "TIM-OS Setup Script"
echo "------------------------------------------------"

# 1. Build C++ Engine
echo "Step 1: Compiling TIM-OS engine..."
make clean && make
echo "C++ compilation successful."
echo

# 2. Setup Python Virtual Environment
echo "Step 2: Setting up Python environment for visualization..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo "Created virtual environment 'venv'."
fi

echo "Installing visualization dependencies..."
./venv/bin/pip install --quiet --upgrade pip
./venv/bin/pip install --quiet numpy pyvista imageio imageio-ffmpeg scipy

echo "Python environment setup successful."
echo

# 3. Usage Guide
echo "------------------------------------------------"
echo "Setup Complete!"
echo "------------------------------------------------"
echo "To run a simulation:"
echo "  ./source/timos -p <opt> -f <mesh> -s <source> -m is -o result.dat"
echo
echo "To open the Interactive Visual Tuner:"
echo "  ./venv/bin/python3 tools/viz_timos.py --mesh <mesh> --data result.dat --grid result.dat.grid --log --interactive --animate"
echo
echo "Example (Half-Sphere Pulse):"
echo "  ./source/timos -f example/half_sphere/spherelens.mesh -s example/half_sphere/sphere_large.source -p example/half_sphere/tissue.opt -m is -o res.dat -T 0.0001 300 -g 3.0 5.0 200 200 -t 8"
echo "------------------------------------------------"
