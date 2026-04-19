from pymol import cmd
from pymol import stored
import sys


def box(dx, dy, dz, cx=None, cy=None, cz=None, pad=None, name="box"):
    """
    Create a box using PyMOL pseudoatoms connected by bonds.
    
    USAGE:
        box dx, dy, dz [, cx, cy, cz, pad [, name]]
    
    ARGUMENTS:
        dx, dy, dz = box dimensions (side lengths)
        cx, cy, cz = box center coordinates (optional)
        pad = distance from the selection to the box walls
        name = name for the box object (default: "box")
    
    NOTES:
        - If an atom is selected and center coordinates are not provided,
          the box will be centered on the selected atom.
        - If no selection and no center coordinates, defaults to origin (0,0,0).
    
    EXAMPLES:
        box 10, 10, 10                    # 10x10x10 box at selected atom or origin
        box 20, 15, 10, 5, 5, 5          # 20x15x10 box centered at (5,5,5)
        box 10, 10, 10, name=mybox       # Custom name for the box object
    """
    
    # Convert dimensions to float
    dx, dy, dz = float(dx), float(dy), float(dz)
    if pad != None:
      pad = float(pad)
    
    # Determine box center
    if cx is None or cy is None or cz is None:
        # Try to get coordinates from selection
        selection = cmd.get_object_list('sele')
        if selection:
            # Get center of mass of selection
            model = cmd.get_model('sele')
            if len(model.atom) > 0:
                cx, cy, cz = 0, 0, 0
                cxmin, cxmax = 1000000,-1000000
                cymin, cymax = 1000000,-1000000
                czmin, czmax = 1000000,-1000000
                for atom in model.atom:
                    cxmin = atom.coord[0] if atom.coord[0] < cxmin else cxmin
                    cymin = atom.coord[1] if atom.coord[1] < cymin else cymin
                    czmin = atom.coord[2] if atom.coord[2] < czmin else czmin
                    cxmax = atom.coord[0] if atom.coord[0] > cxmax else cxmax
                    cymax = atom.coord[1] if atom.coord[1] > cymax else cymax
                    czmax = atom.coord[2] if atom.coord[2] > czmax else czmax
                    cx += atom.coord[0]
                    cy += atom.coord[1]
                    cz += atom.coord[2]
                cx = cx / len(model.atom)
                cy = cy / len(model.atom)
                cz = cz / len(model.atom)
                if pad != None :
                    dx = cxmax - cxmin + 2*pad
                    dy = cymax - cymin + 2*pad
                    dz = czmax - czmin + 2*pad
                #cx = sum(atom.coord[0] for atom in model.atom) / len(model.atom)
                #cy = sum(atom.coord[1] for atom in model.atom) / len(model.atom)
                #cz = sum(atom.coord[2] for atom in model.atom) / len(model.atom)
                print(f"Box centered on selection: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
            else:
                cx, cy, cz = 0.0, 0.0, 0.0
                print("No atoms in selection. Using origin (0, 0, 0)")
        else:
            view = cmd.get_view()
            cx, cy, cz = view[12], view[13], view[14]
            print("No selection. Using origin (0, 0, 0)")
    else:
        cx, cy, cz = float(cx), float(cy), float(cz)
        print(f"Box centered at: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    
    # Calculate half dimensions
    hx, hy, hz = dx / 2.0, dy / 2.0, dz / 2.0
    
    # Calculate the 8 corner vertices
    vertices = [
        (cx - hx, cy - hy, cz - hz),  # 0: ---
        (cx + hx, cy - hy, cz - hz),  # 1: +--
        (cx + hx, cy + hy, cz - hz),  # 2: ++-
        (cx - hx, cy + hy, cz - hz),  # 3: -+-
        (cx - hx, cy - hy, cz + hz),  # 4: --+
        (cx + hx, cy - hy, cz + hz),  # 5: +-+
        (cx + hx, cy + hy, cz + hz),  # 6: +++
        (cx - hx, cy + hy, cz + hz),  # 7: -++
    ]
    
    # Define edges (pairs of vertex indices to connect)
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # bottom face
        (4, 5), (5, 6), (6, 7), (7, 4),  # top face
        (0, 4), (1, 5), (2, 6), (3, 7),  # vertical edges
    ]
    
    # Create pseudoatoms for vertices
    cmd.delete(name)  # Remove if exists
    for i, (x, y, z) in enumerate(vertices):
        cmd.pseudoatom(name, pos=[x, y, z], name=f"V{i}")
    
    # Create bonds between vertices
    for v1, v2 in edges:
        cmd.bond(f"{name} and name V{v1}", f"{name} and name V{v2}")
    
    # Set display properties
    cmd.show("sticks", name)
    cmd.color("white", name)
    cmd.set("stick_radius", 0.1, name)
    
    # Calculate and print box information
    x_min, x_max = cx - hx, cx + hx
    y_min, y_max = cy - hy, cy + hy
    z_min, z_max = cz - hz, cz + hz
    
    print("\n" + "=" * 50)
    print("BOX INFORMATION")
    print("=" * 50)
    print(f"Center:      ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"Dimensions:  {dx:.3f} x {dy:.3f} x {dz:.3f}")
    print(f"X extent:    {x_min:.3f} to {x_max:.3f}")
    print(f"Y extent:    {y_min:.3f} to {y_max:.3f}")
    print(f"Z extent:    {z_min:.3f} to {z_max:.3f}")
    print("=" * 50 + "\n")
    
    return name


# Extend PyMOL command interface
cmd.extend("box", box)

print("Box command loaded successfully!")
print("Usage: box dx, dy, dz [, cx, cy, cz [, name]]")
