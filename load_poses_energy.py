from pymol import cmd


def load_poses_energy(pdbqt_file):
    """
    Load Vina docking poses and display their binding energies as labels.

    Usage in PyMOL:
        load_poses_energy P17_docked_e20_1m52.pdbqt
        load_poses_energy imatinib_docked_e20.pdbqt
    """
    with open(pdbqt_file, "r") as fd:
        lines = fd.read().splitlines()

    energies = []
    for line in lines:
        if line[0:18] == 'REMARK VINA RESULT':
            energies.append(float(line[20:32]))

    cmd.set('label_digits', 4)
    cmd.set('label_relative_mode', 2)
    cmd.load(pdbqt_file, 'obj')
    num_states = cmd.count_states('obj')

    for i in range(num_states, 0, -1):
        patom = '_ene_' + str(i)
        cmd.pseudoatom(patom, 'obj', state=1)
        cmd.set('label_screen_point', [100, 900 + (10 - i) * 50, 0], patom)
        cmd.label(patom, str(energies[i - 1]))
        cmd.set_title('obj', i, energies[i - 1])
    cmd.hide('nonbonded', '_ene_*')

    print(f"Loaded {num_states} poses from {pdbqt_file}")


cmd.extend("load_poses_energy", load_poses_energy)