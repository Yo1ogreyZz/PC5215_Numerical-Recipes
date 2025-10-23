#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import sys

def main():
    filepath = 'Assignment3.ipynb'

    # Read notebook
    with open(filepath, 'r', encoding='utf-8') as f:
        nb = json.load(f)

    # Experimental data to update
    updates = {
        'R_eq': 1.5467,
        'R_eq_ang': 0.8185,
        'E_min': -1.121295,
        'E_min_err': 0.003327,
        'E_bind': 3.3006,
        'E_bind_err': 0.0924,
        'E_eq': -1.108519,
        'E_eq_err': 0.003337,
        'E_eq_ev': -30.1643,
        'E_eq_err_ev': 0.0908,
        'acc_eq': 0.5351,
        'sample_std': 0.425767,
    }

    cells_found = 0

    # Iterate through cells to find and update print statements
    for idx, cell in enumerate(nb['cells']):
        if cell['cell_type'] != 'code':
            continue

        source = ''.join(cell['source'])

        # Pattern 1: Find the cell with R_eq, E_min, E_bind prints
        if 'print(f"R_eq:' in source and 'R_eq_ang' in source and 'E_bind_ev' in source:
            print(f"Found Pattern 1 in cell {idx}")
            cell['source'] = f'''print(f"{{'R_eq:':<25}}{{R_eq:8.4f}} Bohr ({{R_eq_ang:8.4f}} Å)")
print(f"{{'E_min:':<25}}{{E_min:10.6f}} ± {{E_min_err:0.6f}} Hartree")
print(f"{{'E_bind:':<25}}{{E_bind_ev:8.4f}} ± {{E_bind_err_ev:0.4f}} eV")'''
            cells_found += 1

        # Pattern 2: Find the cell with Energy, Acceptance rate, Sample std dev prints
        if 'print(f"Energy:' in source and 'E_eq*HARTREE_TOEV' in source:
            print(f"Found Pattern 2 in cell {idx}")
            cell['source'] = f'''print(f"{{'Energy:':<25}}{{E_eq:10.6f}} ± {{E_eq_err:0.6f}} Hartree ({{E_eq*HARTREE_TOEV:8.4f}} ± {{E_eq_err*HARTREE_TOEV:0.4f}} eV)")
print(f"{{'Acceptance rate:':<25}}{{acc_eq*100:6.2f}}%")
print(f"{{'Sample std dev:':<25}}{{np.std(energies_eq):0.6f}} Hartree")'''
            cells_found += 1

    # Write back
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)

    print(f"\nNotebook updated! Found and updated {cells_found} cells.")
    print("\nYour experimental data:")
    for key, val in updates.items():
        print(f"  {key}: {val}")

if __name__ == '__main__':
    main()

