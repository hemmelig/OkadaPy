"""

"""

import numpy as np
import pandas as pd

from okada import evaluate, read_model


model = read_model("conorbacon_ds01.inp")

# Compute the OkadaPy displacement
displacement = evaluate(model, 0.0, "displacement", threads=24).displacement

coulomb3402_disp_benchmark = pd.read_csv(
    "benchmarks/conorbacon_ds01_displacement.cou",
    sep="\s+",
    skiprows=3,
    header=None,
)

print("\nComparing computed displacements...")
print(
    "   ...x-displacement matching: "
    f"{np.allclose(displacement[:, :, 0].flatten(), coulomb3402_disp_benchmark[3].values)}"
)
print(
    "   ...y-displacement matching: "
    f"{np.allclose(displacement[:, :, 1].flatten(), coulomb3402_disp_benchmark[4].values)}"
)
print(
    "   ...z-displacement matching: "
    f"{np.allclose(displacement[:, :, 2].flatten(), coulomb3402_disp_benchmark[5].values)}"
)
print("...complete.\n")

# Compute the OkadaPy strain
strain = evaluate(model, 0.0, "strain", threads=24).strain

coulomb3402_strain_benchmark = pd.read_csv(
    "benchmarks/conorbacon_ds01_strain.cou",
    sep="\s+",
    skiprows=3,
    header=None,
)

print("\nComparing computed strains...")
print(
    "   ...exx matching: "
    f"{np.allclose(strain[:, :, 0].flatten(), coulomb3402_strain_benchmark[3].values)}"
)
print(
    "   ...eyy matching: "
    f"{np.allclose(strain[:, :, 1].flatten(), coulomb3402_strain_benchmark[4].values)}"
)
print(
    "   ...ezz matching: "
    f"{np.allclose(strain[:, :, 2].flatten(), coulomb3402_strain_benchmark[5].values)}"
)
print(
    "   ...exy matching: "
    f"{np.allclose(strain[:, :, 5].flatten(), coulomb3402_strain_benchmark[8].values)}"
)
print(
    "   ...exz matching: "
    f"{np.allclose(strain[:, :, 4].flatten(), coulomb3402_strain_benchmark[7].values)}"
)
print(
    "   ...eyz matching: "
    f"{np.allclose(strain[:, :, 3].flatten(), coulomb3402_strain_benchmark[6].values)}"
)
print("...complete.\n")

# Compute the OkadaPy stress
stress = evaluate(model, 0.0, "stress", threads=24).stress

coulomb3402_stress_benchmark = pd.read_csv(
    "benchmarks/conorbacon_ds01_stress.cou",
    sep="\s+",
    skiprows=3,
    header=None,
)

print("\nComparing computed stresses...")
print(
    "   ...sxx matching: "
    f"{np.allclose(stress[:, :, 0].flatten(), coulomb3402_stress_benchmark[3].values)}"
)
print(
    "   ...syy matching: "
    f"{np.allclose(stress[:, :, 1].flatten(), coulomb3402_stress_benchmark[4].values, rtol=1e-04)}"
)
print(
    "   ...szz matching: "
    f"{np.allclose(stress[:, :, 2].flatten(), coulomb3402_stress_benchmark[5].values)}"
)
print(
    "   ...sxy matching: "
    f"{np.allclose(stress[:, :, 5].flatten(), coulomb3402_stress_benchmark[8].values)}"
)
print(
    "   ...sxz matching: "
    f"{np.allclose(stress[:, :, 4].flatten(), coulomb3402_stress_benchmark[7].values)}"
)
print(
    "   ...syz matching: "
    f"{np.allclose(stress[:, :, 3].flatten(), coulomb3402_stress_benchmark[6].values)}"
)
print("...complete.\n")
