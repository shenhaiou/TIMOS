import sys
import numpy as np
import pyvista as pv
import argparse
import os

def parse_mesh(mesh_path):
    print(f"Reading mesh: {mesh_path}")
    with open(mesh_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    num_nodes = int(lines[0])
    num_elems = int(lines[1])
    nodes = np.array([[float(x) for x in lines[i].split()] for i in range(2, 2 + num_nodes)])
    cells = []
    cell_types = []
    for i in range(2 + num_nodes, 2 + num_nodes + num_elems):
        parts = [int(x) for x in lines[i].split()]
        cells.append([4, parts[0]-1, parts[1]-1, parts[2]-1, parts[3]-1])
        cell_types.append(pv.CellType.TETRA)
    return pv.UnstructuredGrid(np.hstack(cells), cell_types, nodes)

def parse_dat(dat_path, num_elems, num_boundary_trigs):
    print(f"Reading data: {dat_path}")
    fluence_int, fluence_surf, num_steps = None, None, 1
    with open(dat_path, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('1 '):
            parts = line.split()
            count, steps = int(parts[1]), int(parts[2])
            if count == num_boundary_trigs:
                num_steps = steps
                fluence_surf = np.array([[float(x) for x in lines[i+1+j].split()[4:]] for j in range(count)])
        elif line.startswith('2 '):
            parts = line.split()
            count, steps = int(parts[1]), int(parts[2])
            if count == num_elems:
                num_steps = steps
                fluence_int = np.array([[float(x) for x in lines[i+1+j].split()[5:]] for j in range(count)])
    return fluence_int, fluence_surf, num_steps

def parse_grid(grid_path):
    print(f"Reading grid data: {grid_path}")
    with open(grid_path, 'r') as f:
        lines = f.readlines()
    nr, ny, nt = 0, 0, 0
    start_idx = -1
    for i, line in enumerate(lines):
        if line.startswith('3 '):
            parts = line.split()
            nr, ny, nt = int(parts[1]), int(parts[2]), int(parts[3])
            start_idx = i + 1
            break
    data = np.zeros((nr, ny, nt))
    for i in range(nr * ny):
        parts = lines[start_idx + i].split()
        data[int(parts[0]), int(parts[1]), :] = [float(x) for x in parts[2:]]
    return data, nr, ny, nt

def main():
    parser = argparse.ArgumentParser(description="TIM-OS Visualization Tool")
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--grid", help="Path to .grid result file")
    parser.add_argument("--surface", action="store_true")
    parser.add_argument("--slice-z", type=float)
    parser.add_argument("--slice-y", type=float)
    parser.add_argument("--slice-x", type=float)
    parser.add_argument("--output")
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--animate", action="store_true")
    args = parser.parse_args()

    grid_mesh = parse_mesh(args.mesh)
    num_bt = 0
    with open(args.data, 'r') as f:
        for line in f:
            if line.startswith('1 '): num_bt = int(line.split()[1]); break

    fluence_int, fluence_surf, num_steps = parse_dat(args.data, grid_mesh.n_cells, num_bt)
    
    # OLD BLUE SCHEMA: Jet-like (Blue to Red)
    cmap = "jet" 

    if args.grid:
        data, nr, ny, nt = parse_grid(args.grid)
        num_steps = nt
        rmax, ymax = 3.0, 5.0
        with open(args.grid, 'r') as f:
            for l in f.readlines()[:10]:
                if "Rmax=" in l:
                    parts = l.replace("%","").strip().split()
                    for p in parts:
                        if "Rmax=" in p: rmax = float(p.split('=')[1])
                        if "Ymax=" in p: ymax = float(p.split('=')[1])
                    break
        mesh = pv.ImageData(dimensions=(2*nr, ny, 1), spacing=(rmax/nr, ymax/ny, 1))
        mesh.origin = (-rmax, 0, 0)
        def get_step_data(step):
            d = np.zeros((2*nr, ny))
            d[:nr, :] = np.flipud(data[:, :, step])
            d[nr:, :] = data[:, :, step]
            return np.log10(d.T + 1e-18) if args.log else d.T
    else:
        mesh = grid_mesh
        fluence = fluence_surf if args.surface else fluence_int
        def get_step_data(step):
            d = fluence[:, step]
            return np.log10(d + 1e-18) if args.log else d

    v_max = np.nanmax(get_step_data(num_steps//4 if num_steps>1 else 0))
    v_min = v_max - 8 if args.log else 0

    plotter = pv.Plotter(off_screen=True if args.output else False)
    d0 = get_step_data(0)
    if args.grid: mesh.point_data["Fluence"] = d0.flatten()
    else: mesh.cell_data["Fluence"] = d0

    display_mesh = mesh
    if not args.grid and not args.surface:
        has_slice = any(v is not None for v in [args.slice_x, args.slice_y, args.slice_z])
        if has_slice:
            if args.slice_z is not None: display_mesh = mesh.slice(normal='z', origin=(0, 0, args.slice_z))
            elif args.slice_y is not None: display_mesh = mesh.slice(normal='y', origin=(0, args.slice_y, 0))
            else: display_mesh = mesh.slice(normal='x', origin=(args.slice_x, 0, 0))
        else: display_mesh = mesh.threshold(1e-15 if not args.log else -17, scalars="Fluence")

    plotter.add_mesh(display_mesh, scalars="Fluence", cmap=cmap, clim=[v_min, v_max], scalar_bar_args={"title": "Log10(Fluence)" if args.log else "Fluence"})
    plotter.view_xy() # Flat view
    
    if args.animate and num_steps > 1:
        video_path = args.output if args.output and args.output.endswith('.mp4') else "animation.mp4"
        plotter.open_movie(video_path)
        print(f"Generating animation to {video_path}...")
        for i in range(num_steps):
            step_data = get_step_data(i)
            if args.grid: display_mesh.point_data["Fluence"] = step_data.flatten()
            else: display_mesh.cell_data["Fluence"] = step_data
            plotter.write_frame()
        plotter.close()
    elif args.output:
        plotter.screenshot(args.output)
    else:
        plotter.show()

if __name__ == "__main__":
    main()
