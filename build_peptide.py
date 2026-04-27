from pymol import cmd

def build_peptide_uniform(sequence, name="peptide", 
                          phi=-60.0, psi=-45.0, omega=180.0):
    """
    Builds a peptide with the SAME phi, psi, and omega angles repeated for every residue.
    
    Parameters:
    sequence : str     One-letter amino acid sequence (e.g. "AAAAAAAAAA")
    name     : str     Name of the object in PyMOL
    phi      : float   Phi angle in degrees (default -60°)
    psi      : float   Psi angle in degrees (default -45°)
    omega    : float   Omega (peptide bond) angle in degrees (default 180° trans)
    """
    if len(sequence) < 1:
        print("Error: Sequence must have at least one residue.")
        return
    
    # Build an extended peptide first
    cmd.fab(sequence, name, ss=0)
    
    # Apply the same phi/psi/omega to all residues
    for i in range(len(sequence)):
        resi = i + 1
        
        # Set Phi (except for the first residue)
        if i > 0:
            cmd.set_dihedral(f"{name}///{resi-1}/C", 
                             f"{name}///{resi}/N", 
                             f"{name}///{resi}/CA", 
                             f"{name}///{resi}/C", 
                             phi)
        
        # Set Psi (except for the last residue)
        if i < len(sequence) - 1:
            cmd.set_dihedral(f"{name}///{resi}/N", 
                             f"{name}///{resi}/CA", 
                             f"{name}///{resi}/C", 
                             f"{name}///{resi+1}/N", 
                             psi)
        
        # Set Omega (except for the last residue)
        if i < len(sequence) - 1:
            cmd.set_dihedral(f"{name}///{resi}/CA", 
                             f"{name}///{resi}/C", 
                             f"{name}///{resi+1}/N", 
                             f"{name}///{resi+1}/CA", 
                             omega)
    
    # Finalize display
    cmd.center(name)
    cmd.show("sticks", name)
    cmd.show("cartoon", name)
    cmd.util.cbc(name)          # color by chain
    cmd.zoom(name, 1.2)
    
    print(f"Peptide '{name}' built successfully!")
    print(f"   Length     : {len(sequence)} residues")
    print(f"   Phi (φ)    : {phi}° (repeated)")
    print(f"   Psi (ψ)    : {psi}° (repeated)")
    print(f"   Omega (ω)  : {omega}° (repeated)")


# Register the command so you can call it directly in PyMOL
cmd.extend("build_peptide_uniform", build_peptide_uniform)