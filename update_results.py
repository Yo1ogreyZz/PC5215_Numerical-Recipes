import json

# Read the notebook
with open('Assignment3.ipynb', encoding='utf-8') as f:
    nb = json.load(f)

# Your experimental data
R_eq = 1.5467
R_eq_ang = 0.8185
E_min = -1.121295
E_min_err = 0.003327
E_bind = 3.3006
E_bind_err = 0.0924
E_eq = -1.108519
E_eq_err = 0.003337
E_eq_ev = -30.1643
E_eq_err_ev = 0.0908
acc_eq = 0.5351
sample_std = 0.425767

# Find and update relevant cells
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])

        # Update the first print statement (R_eq, E_min, E_bind)
        if 'R_eq: {R_eq:.4f}' in source and 'E_bind_ev' in source:
            new_source = f'''print(f"{{'R_eq:':<25}}{{R_eq:8.4f}} Bohr ({{R_eq_ang:8.4f}} Å)")
print(f"{{'E_min:':<25}}{{E_min:10.6f}} ± {{E_min_err:0.6f}} Hartree")
print(f"{{'E_bind:':<25}}{{E_bind_ev:8.4f}} ± {{E_bind_err_ev:0.4f}} eV")'''
            cell['source'] = new_source

        # Update the equilibrium detailed analysis print
        if 'Energy: {E_eq:.6f}' in source and 'Acceptance rate' in source:
            new_source = f'''print(f"{{'Energy:':<25}}{{E_eq:10.6f}} ± {{E_eq_err:0.6f}} Hartree ({{E_eq_ev:8.4f}} ± {{E_eq_err_ev:0.4f}} eV)")
print(f"{{'Acceptance rate:':<25}}{{acc_eq*100:6.2f}}%")
print(f"{{'Sample std dev:':<25}}{{sample_std:0.6f}} Hartree")'''
            cell['source'] = new_source

# Write the notebook back
with open('Assignment3.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print("Notebook updated successfully!")
print(f"\nUpdated values:")
print(f"R_eq: {R_eq:.4f} Bohr ({R_eq_ang:.4f} Å)")
print(f"E_min: {E_min:.6f} ± {E_min_err:.6f} Hartree")
print(f"E_bind: {E_bind:.4f} ± {E_bind_err:.4f} eV")
print(f"Energy: {E_eq:.6f} ± {E_eq_err:.6f} Hartree ({E_eq_ev:.4f} ± {E_eq_err_ev:.4f} eV)")
print(f"Acceptance rate: {acc_eq*100:.2f}%")
print(f"Sample std dev: {sample_std:.6f} Hartree")

