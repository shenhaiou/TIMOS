import sys
import numpy as np
import pyvista as pv
import argparse
import os

def parse_mesh(mesh_path):
    with open(mesh_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    num_nodes, num_elems = int(lines[0]), int(lines[1])
    nodes = np.array([[float(x) for x in lines[i].split()] for i in range(2, 2 + num_nodes)])
    cells = []
    for i in range(2 + num_nodes, 2 + num_nodes + num_elems):
        parts = [int(x) for x in lines[i].split()]
        cells.append([4, parts[0]-1, parts[1]-1, parts[2]-1, parts[3]-1])
    return pv.UnstructuredGrid(np.hstack(cells), [pv.CellType.TETRA]*num_elems, nodes)

def parse_grid(grid_path):
    with open(grid_path, 'r') as f:
        lines = f.readlines()
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--grid")
    parser.add_argument("--output")
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--animate", action="store_true")
    parser.add_argument("--smooth", type=float, default=0.0)
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--drange", type=float, help="Dynamic range in log units (e.g. 4.0)")
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--preview-frame", type=int, help="Frame to show in static/interactive preview")
    parser.add_argument("--interactive", action="store_true", help="Tune color scale manually")
    args = parser.parse_args()

    data, nr, ny, nt = parse_grid(args.grid)
    if args.max_steps and args.max_steps < nt: nt = args.max_steps
    if args.smooth > 0:
        from scipy.ndimage import gaussian_filter
        for t in range(nt): data[:, :, t] = gaussian_filter(data[:, :, t], sigma=args.smooth)

    rmax, ymax = 3.0, 5.0
    with open(args.grid, 'r') as f:
        for l in f.readlines()[:10]:
            if "Rmax=" in l:
                parts = l.replace("%","").strip().split()
                for p in parts:
                    if "Rmax=" in p: rmax = float(p.split('=')[1])
                    if "Ymax=" in p: ymax = float(p.split('=')[1])

    mesh = pv.ImageData(dimensions=(2*nr, ny, 1), spacing=(rmax/nr, ymax/ny, 1))
    mesh.origin = (-rmax, 0, 0)
    
    def get_step_data(t):
        d = np.zeros((2*nr, ny))
        d[:nr, :] = np.flipud(data[:, :, t]); d[nr:, :] = data[:, :, t]
        return np.log10(d.T + 1e-18) if args.log else d.T

    # Find global max
    actual_max = -1e20
    for t in range(nt):
        m = np.nanmax(get_step_data(t))
        if m > actual_max: actual_max = m
    
    v_max = args.vmax if args.vmax is not None else actual_max
    v_min = args.vmin if args.vmin is not None else (v_max - (args.drange if args.drange else 4.0))

    if args.interactive:
        print(f"\n--- Interactive Color Tuning ---")
        print(f"Data Peak: {actual_max:.2f}")
        while True:
            print(f"Current Range: {v_min:.2f} to {v_max:.2f}")
            plotter = pv.Plotter(title="TIM-OS Tuning Preview (Close to accept)")
            preview_idx = args.preview_frame if args.preview_frame is not None else nt // 4
            mesh.point_data["Fluence"] = get_step_data(preview_idx).flatten()
            plotter.add_mesh(mesh, scalars="Fluence", cmap="jet", clim=[v_min, v_max])
            plotter.view_xy()
            plotter.show()
            
            resp = input("Accept range? (y/n) or enter new 'vmin vmax' or 'drange': ")
            if resp.lower() == 'y' or resp == '': break
            try:
                parts = resp.split()
                if len(parts) == 1:
                    args.drange = float(parts[0])
                    v_min = v_max - args.drange
                else:
                    v_min, v_max = float(parts[0]), float(parts[1])
            except:
                print("Invalid input. Try again.")

    print(f"Final Range: {v_min} to {v_max}")

    plotter = pv.Plotter(off_screen=True if args.output else False)
    # Use preview_frame for static image if provided
    d0 = get_step_data(args.preview_frame if args.preview_frame is not None else 0)
    mesh.point_data["Fluence"] = d0.flatten()
    actor = plotter.add_mesh(mesh, scalars="Fluence", cmap="jet", clim=[v_min, v_max])
    plotter.view_xy()
    
    if args.animate and nt > 1:
        video_path = args.output if args.output else "animation.mp4"
        plotter.open_movie(video_path)
        for i in range(nt):
            mesh.point_data["Fluence"] = get_step_data(i).flatten()
            plotter.write_frame()
        plotter.close()
    elif args.output: plotter.screenshot(args.output)
    else: plotter.show()

if __name__ == "__main__": main()
