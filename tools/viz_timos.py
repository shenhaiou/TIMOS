import sys
import numpy as np
import pyvista as pv
import argparse
import os

def parse_mesh(mesh_path):
    print(f"Reading mesh: {mesh_path}")
    with open(mesh_path, 'r') as f:
        # Skip comments and empty lines
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    
    num_nodes = int(lines[0])
    num_elems = int(lines[1])
    
    # Parse nodes
    nodes = []
    for i in range(2, 2 + num_nodes):
        nodes.append([float(x) for x in lines[i].split()])
    nodes = np.array(nodes)
    
    # Parse elements
    cells = []
    cell_types = []
    for i in range(2 + num_nodes, 2 + num_nodes + num_elems):
        parts = [int(x) for x in lines[i].split()]
        # parts: n1, n2, n3, n4, medIdx
        # VTK format for UnstructuredGrid: [n_points, p1, p2, p3, p4]
        # Mesh indices are 1-based, convert to 0-based
        cells.append([4, parts[0]-1, parts[1]-1, parts[2]-1, parts[3]-1])
        cell_types.append(pv.CellType.TETRA)
    
    cells = np.hstack(cells)
    grid = pv.UnstructuredGrid(cells, cell_types, nodes)
    return grid

def parse_dat(dat_path, num_elems, num_boundary_trigs):
    print(f"Reading data: {dat_path}")
    fluence_int = None
    fluence_surf = None
    num_steps = 1
    
    with open(dat_path, 'r') as f:
        lines = f.readlines()
    
    # Locate sections
    for i, line in enumerate(lines):
        if line.startswith('1 '): # Surface
            parts = line.split()
            count = int(parts[1])
            num_steps = int(parts[2])
            if count == num_boundary_trigs:
                fluence_surf = np.zeros((count, num_steps))
                for j in range(count):
                    data = lines[i + 1 + j].split()
                    # n1 n2 n3 area flu1 flu2 ...
                    fluence_surf[j, :] = [float(x) for x in data[4:]]
        
        elif line.startswith('2 '): # Internal
            parts = line.split()
            count = int(parts[1])
            num_steps = int(parts[2])
            if count == num_elems:
                fluence_int = np.zeros((count, num_steps))
                for j in range(count):
                    data = lines[i + 1 + j].split()
                    # n1 n2 n3 n4 vol flu1 flu2 ...
                    fluence_int[j, :] = [float(x) for x in data[5:]]
                    
    return fluence_int, fluence_surf, num_steps

def get_surface_mesh(mesh_path, dat_path):
    # We need nodes from mesh, but triangles from .dat since .mesh doesn't list them easily
    with open(mesh_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    num_nodes = int(lines[0])
    nodes = np.array([[float(x) for x in lines[i].split()] for i in range(2, 2 + num_nodes)])
    
    with open(dat_path, 'r') as f:
        lines = f.readlines()
    
    faces = []
    for i, line in enumerate(lines):
        if line.startswith('1 '):
            count = int(line.split()[1])
            for j in range(count):
                parts = lines[i + 1 + j].split()
                # n1 n2 n3 (1-based)
                faces.append([3, int(parts[0])-1, int(parts[1])-1, int(parts[2])-1])
            break
    
    if not faces:
        return None
    
    return pv.PolyData(nodes, np.hstack(faces))

def main():
    parser = argparse.ArgumentParser(description="TIM-OS Visualization Tool")
    parser.add_argument("--mesh", required=True, help="Path to .mesh file")
    parser.add_argument("--data", required=True, help="Path to .dat result file")
    parser.add_argument("--surface", action="store_true", help="Visualize surface fluence")
    parser.add_argument("--slice-z", type=float, help="Coordinate to slice at Z")
    parser.add_argument("--slice-y", type=float, help="Coordinate to slice at Y")
    parser.add_argument("--slice-x", type=float, help="Coordinate to slice at X")
    parser.add_argument("--output", help="Save plot to image file (or .mp4 for animations)")
    parser.add_argument("--log", action="store_true", help="Use log scale for fluence")
    parser.add_argument("--animate", action="store_true", help="Create animation for time-domain data")
    
    args = parser.parse_args()

    if not os.path.exists(args.mesh) or not os.path.exists(args.data):
        print("Error: Mesh or Data file not found.")
        return

    # Initial parse to get counts
    grid = parse_mesh(args.mesh)
    
    # Get boundary trig count from .dat
    num_boundary_trigs = 0
    with open(args.data, 'r') as f:
        for line in f:
            if line.startswith('1 '):
                num_boundary_trigs = int(line.split()[1])
                break

    fluence_int, fluence_surf, num_steps = parse_dat(args.data, grid.n_cells, num_boundary_trigs)
    
    if args.surface:
        mesh = get_surface_mesh(args.mesh, args.data)
        fluence = fluence_surf
        if mesh is None:
            print("Error: No surface data found in .dat file.")
            return
    else:
        mesh = grid
        fluence = fluence_int

    if fluence is None:
        print("Error: Requested data type not found in .dat file.")
        return

    print(f"Mesh bounds: {mesh.bounds}")
    print(f"Time steps: {num_steps}")

    def get_step_data(step):
        d = fluence[:, step]
        if args.log:
            return np.log10(d + 1e-15)
        return d

    label = "Log10(Fluence)" if args.log else "Fluence"

    # Check if any slice is requested
    has_slice = any(v is not None for v in [args.slice_x, args.slice_y, args.slice_z])

    if args.animate and num_steps > 1:
        plotter = pv.Plotter(off_screen=True if args.output else False)
        
        # Initial frame
        mesh.cell_data["Fluence"] = get_step_data(0)
        
        if not args.surface and has_slice:
            if args.slice_z is not None:
                display_mesh = mesh.slice(normal='z', origin=(0, 0, args.slice_z))
            elif args.slice_y is not None:
                display_mesh = mesh.slice(normal='y', origin=(0, args.slice_y, 0))
            else:
                display_mesh = mesh.slice(normal='x', origin=(args.slice_x, 0, 0))
        else:
            display_mesh = mesh

        # Use auto-clim based on log or linear data
        v_min = np.nanmin(get_step_data(0))
        v_max = np.nanmax(fluence if not args.log else np.log10(fluence + 1e-15))
        
        actor = plotter.add_mesh(display_mesh, scalars="Fluence", cmap="viridis", 
                         scalar_bar_args={"title": label}, clim=[v_min, v_max])
        plotter.add_axes()

        video_path = args.output if args.output and args.output.endswith('.mp4') else "animation.mp4"
        plotter.open_movie(video_path)
        
        print(f"Generating animation to {video_path}...")
        for i in range(num_steps):
            mesh.cell_data["Fluence"] = get_step_data(i)
            
            if not args.surface and has_slice:
                 if args.slice_z is not None:
                    new_slice = mesh.slice(normal='z', origin=(0, 0, args.slice_z))
                 elif args.slice_y is not None:
                    new_slice = mesh.slice(normal='y', origin=(0, args.slice_y, 0))
                 else:
                    new_slice = mesh.slice(normal='x', origin=(args.slice_x, 0, 0))
                 display_mesh.shallow_copy(new_slice)
            
            # Since we used display_mesh in add_mesh, shallow_copy updates the plotter mesh
            plotter.write_frame()
        plotter.close()
        print("Done.")
    else:
        # Standard Plotting
        plotter = pv.Plotter(off_screen=True if args.output else False)
        mesh.cell_data["Fluence"] = get_step_data(0)

        if args.surface:
            display_mesh = mesh
        elif args.slice_z is not None:
            display_mesh = mesh.slice(normal='z', origin=(0, 0, args.slice_z))
        elif args.slice_y is not None:
            display_mesh = mesh.slice(normal='y', origin=(0, args.slice_y, 0))
        elif args.slice_x is not None:
            display_mesh = mesh.slice(normal='x', origin=(args.slice_x, 0, 0))
        else:
            display_mesh = mesh.threshold(1e-12 if not args.log else -14, scalars="Fluence")

        if display_mesh.n_points == 0:
            print("Warning: Resulting mesh is empty. Showing full mesh.")
            display_mesh = mesh

        plotter.add_mesh(display_mesh, scalars="Fluence", cmap="viridis", 
                         scalar_bar_args={"title": label})
        plotter.add_axes()
        
if __name__ == "__main__":
    main()
