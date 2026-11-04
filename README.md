# CR3BP Orbital Family Continuation & 3D Visualization

A Python framework for generating and tracing orbital families in the Circular Restricted Three-Body Problem (CR3BP) using **Pseudo-Arclength Continuation (PAC)** powered by the [ASSET](https://github.com/ASDL-ASRL/ASSET) library (`asset_asrl`).

---

## 🚀 Features
* **Adaptive Pseudo-Arclength Continuation:** Automatically scales step sizes (`dS`) to navigate high-gradient gravitational fields near primary bodies (Earth and Moon).
* **Robust Convergence Handling:** Built-in retry loops and gentle step-scaling mechanisms to prevent crashes in highly nonlinear regions.
* **Surface Impact Detection:** Real-time 3D distance tracking against planetary radii to automatically halt continuation upon reaching user-defined surface constraints.
* **Universal 3D Plotting:** Renders 3D trajectory families with Earth, Moon, and collinear Lagrange points (`L1`–`L3`).

---

## 📂 Project Structure
```text
├── PAC_Utilities.py      # Core FamilySolver class (Shooting & PAC continuation)
├── l1_halo_generation.py # Generation script for L1/L2 halo families
├── l3_planar_generation.py # Generation script for L3 planar families
├── Plotter.py            # 3D visualization tool for saved orbit pickle files
└── README.md


📦 Prerequisites & Environment

Ensure your Python environment (such as an Anaconda environment) has the following packages installed:

    asset_asrl

    numpy

    matplotlib

    scipy

⚙️ How to Run
1. Generate an Orbital Family

Run a generation script to compute the orbit family and export the results to a .pkl file:
Bash

python l3_planar_generation.py

(This will generate and save your trajectory data, e.g., l3_planar_family.pkl or l1_halo_family.pkl)
2. Visualize the Family

Run the plotter script to render the 3D orbit visualization and check closest approach telemetry:
Bash

python Plotter.py

🛠️ Built With

    ASSET (asset_asrl) — Trajectory optimization and optimal control transcription.

    NumPy — Vectorized math and Euclidean distance calculations.

    Matplotlib — 3D space rendering.
