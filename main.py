# ============================================================
# MONTE CARLO RADIATION TRANSPORT
# KLEIN-NISHINA REAL
# ============================================================

# ============================================================
# IMPORTS
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm

# ============================================================
# PASTAS
# ============================================================

os.makedirs("outputs/images", exist_ok=True)
os.makedirs("outputs/datasets", exist_ok=True)

# ============================================================
# CONSTANTES FÍSICAS
# ============================================================

# Energia de repouso do elétron
# keV

ELECTRON_REST_ENERGY = 511.0

# Raio clássico do elétron

CLASSICAL_ELECTRON_RADIUS = 2.817e-13

# ============================================================
# MATERIAL
# ============================================================

class Material:

    def __init__(

        self,

        name,

        density,

        atomic_number,

        mu_total
    ):

        self.name = name

        self.density = density

        self.atomic_number = atomic_number

        self.mu_total = mu_total


# ============================================================
# MATERIAIS
# ============================================================

MATERIALS = {

    "water": Material(

        "Water",

        density=1.0,

        atomic_number=7,

        mu_total=0.15
    ),

    "tissue": Material(

        "Tissue",

        density=1.05,

        atomic_number=7.4,

        mu_total=0.20
    ),

    "aluminum": Material(

        "Aluminum",

        density=2.7,

        atomic_number=13,

        mu_total=0.35
    ),

    "lead": Material(

        "Lead",

        density=11.34,

        atomic_number=82,

        mu_total=1.20
    )
}

# ============================================================
# PROBABILIDADES FÍSICAS
# ============================================================

def calculate_interaction_probabilities(

    energy,

    material
):

    Z = material.atomic_number

    # --------------------------------------------------------
    # FOTOELÉTRICO
    #
    # Aproximação:
    #
    # ~ Z^4 / E^3
    # --------------------------------------------------------

    photoelectric = (

        (Z ** 4)

        /

        ((energy + 1e-6) ** 3)
    )

    # --------------------------------------------------------
    # COMPTON
    #
    # Dependência mais fraca
    # --------------------------------------------------------

    compton = (

        Z

        /

        (energy + 1e-6)
    )

    total = photoelectric + compton

    P_photoelectric = photoelectric / total

    P_compton = compton / total

    return P_photoelectric, P_compton


# ============================================================
# KLEIN-NISHINA
# ============================================================

def klein_nishina_differential_cs(

    energy,

    theta
):

    # --------------------------------------------------------
    # Razão energética
    # --------------------------------------------------------

    epsilon = energy / ELECTRON_REST_ENERGY

    # --------------------------------------------------------
    # Energia espalhada
    # --------------------------------------------------------

    E_ratio = (

        1 /

        (
            1 +
            epsilon *
            (1 - np.cos(theta))
        )
    )

    # --------------------------------------------------------
    # Klein-Nishina
    # --------------------------------------------------------

    term1 = (E_ratio ** 2)

    term2 = (

        E_ratio +

        (1 / E_ratio)

        -

        (np.sin(theta) ** 2)
    )

    dsigma = (

        (CLASSICAL_ELECTRON_RADIUS ** 2)

        / 2

        *

        term1

        *

        term2
    )

    return dsigma


# ============================================================
# AMOSTRAGEM KLEIN-NISHINA
# ============================================================

def sample_klein_nishina_angle(

    energy,

    n_samples=1000
):

    # --------------------------------------------------------
    # Geração angular
    # --------------------------------------------------------

    angles = np.linspace(

        0,

        np.pi,

        n_samples
    )

    probabilities = np.array([

        klein_nishina_differential_cs(
            energy,
            theta
        )

        for theta in angles
    ])

    # --------------------------------------------------------
    # Normalização
    # --------------------------------------------------------

    probabilities /= probabilities.sum()

    # --------------------------------------------------------
    # Amostragem
    # --------------------------------------------------------

    theta = np.random.choice(

        angles,

        p=probabilities
    )

    return theta


# ============================================================
# FÓTON
# ============================================================

class Photon:

    def __init__(

        self,

        x,

        y,

        energy,

        angle
    ):

        self.x = x
        self.y = y

        self.energy = energy

        self.angle = angle

        self.alive = True

        self.path_x = [x]
        self.path_y = [y]

        self.energy_history = [energy]

    # --------------------------------------------------------

    def move(self, distance):

        dx = distance * np.cos(self.angle)
        dy = distance * np.sin(self.angle)

        self.x += dx
        self.y += dy

        self.path_x.append(self.x)
        self.path_y.append(self.y)

        self.energy_history.append(self.energy)

    # --------------------------------------------------------
    # COMPTON COM KLEIN-NISHINA
    # --------------------------------------------------------

    def compton_scatter(self):

        # --------------------------------------------
        # Ângulo físico real
        # --------------------------------------------

        theta = sample_klein_nishina_angle(

            self.energy
        )

        E_initial = self.energy

        # --------------------------------------------
        # Fórmula de Compton
        # --------------------------------------------

        E_scattered = (

            E_initial

            /

            (
                1 +

                (
                    E_initial /
                    ELECTRON_REST_ENERGY
                )

                *

                (
                    1 - np.cos(theta)
                )
            )
        )

        deposited_energy = (

            E_initial - E_scattered
        )

        self.energy = E_scattered

        self.angle += theta

        # --------------------------------------------
        # Critério de parada
        # --------------------------------------------

        if self.energy < 1:

            self.alive = False

        return deposited_energy, theta

    # --------------------------------------------------------
    # FOTOELÉTRICO
    # --------------------------------------------------------

    def photoelectric_absorption(self):

        deposited_energy = self.energy

        self.energy = 0

        self.alive = False

        return deposited_energy


# ============================================================
# MAPA DE DOSE
# ============================================================

class DoseMap:

    def __init__(

        self,

        width,

        height
    ):

        self.width = width
        self.height = height

        self.grid = np.zeros((height, width))

    # --------------------------------------------------------

    def deposit_energy(

        self,

        x,

        y,

        energy
    ):

        ix = int(x)
        iy = int(y)

        if (

            0 <= ix < self.width
            and
            0 <= iy < self.height
        ):

            self.grid[iy, ix] += energy


# ============================================================
# MONTE CARLO
# ============================================================

class MonteCarloSimulation:

    def __init__(

        self,

        material,

        num_particles=30000,

        initial_energy=120.0,

        width=400,

        height=400
    ):

        self.material = material

        self.num_particles = num_particles

        self.initial_energy = initial_energy

        self.width = width
        self.height = height

        self.dose_map = DoseMap(width, height)

        self.dataset = []

        self.trajectories = []

        self.photoelectric_events = 0
        self.compton_events = 0

    # --------------------------------------------------------
    # Caminho livre médio
    # --------------------------------------------------------

    def sample_free_path(self):

        mu = self.material.mu_total

        u = np.random.uniform(1e-10, 1)

        return -np.log(u) / mu

    # --------------------------------------------------------
    # Simular fóton
    # --------------------------------------------------------

    def simulate_photon(self):

        photon = Photon(

            x=self.width // 2,

            y=self.height // 2,

            energy=self.initial_energy,

            angle=np.random.uniform(0, 2*np.pi)
        )

        while photon.alive:

            # --------------------------------------------
            # Movimento
            # --------------------------------------------

            distance = self.sample_free_path()

            photon.move(distance)

            # --------------------------------------------
            # Limites
            # --------------------------------------------

            if (

                photon.x < 0
                or photon.x >= self.width
                or photon.y < 0
                or photon.y >= self.height
            ):

                break

            # --------------------------------------------
            # Probabilidades
            # --------------------------------------------

            P_photoelectric, P_compton = \
                calculate_interaction_probabilities(

                    photon.energy,

                    self.material
                )

            interaction = np.random.random()

            # --------------------------------------------
            # FOTOELÉTRICO
            # --------------------------------------------

            if interaction < P_photoelectric:

                deposited_energy = \
                    photon.photoelectric_absorption()

                self.dose_map.deposit_energy(

                    photon.x,

                    photon.y,

                    deposited_energy
                )

                self.photoelectric_events += 1

                self.dataset.append({

                    "interaction":
                        "photoelectric",

                    "energy_before":
                        deposited_energy,

                    "energy_after":
                        0,

                    "theta":
                        0,

                    "distance":
                        distance,

                    "x":
                        photon.x,

                    "y":
                        photon.y
                })

            # --------------------------------------------
            # COMPTON
            # --------------------------------------------

            else:

                energy_before = photon.energy

                deposited_energy, theta = \
                    photon.compton_scatter()

                self.dose_map.deposit_energy(

                    photon.x,

                    photon.y,

                    deposited_energy
                )

                self.compton_events += 1

                self.dataset.append({

                    "interaction":
                        "compton",

                    "energy_before":
                        energy_before,

                    "energy_after":
                        photon.energy,

                    "theta":
                        np.degrees(theta),

                    "distance":
                        distance,

                    "x":
                        photon.x,

                    "y":
                        photon.y
                })

        self.trajectories.append({

            "x": photon.path_x,

            "y": photon.path_y
        })

    # --------------------------------------------------------

    def run(self):

        for _ in tqdm(range(self.num_particles)):

            self.simulate_photon()

        return {

            "dose_map":
                self.dose_map.grid,

            "dataset":
                pd.DataFrame(self.dataset),

            "trajectories":
                self.trajectories,

            "photoelectric":
                self.photoelectric_events,

            "compton":
                self.compton_events
        }


# ============================================================
# VISUALIZAÇÃO
# ============================================================

class Visualization:

    # --------------------------------------------------------
    # HEATMAP
    # --------------------------------------------------------

    @staticmethod
    def plot_dose_map(dose_map):

        plt.figure(figsize=(12, 10))

        dose_log = np.log1p(dose_map)

        sns.heatmap(

            dose_log,

            cmap="inferno"
        )

        plt.title(

            "Dose Map (Log Scale)",

            fontsize=18
        )

        plt.xlabel("X")
        plt.ylabel("Y")

        plt.savefig(

            "outputs/images/dose_map.png",

            dpi=300,

            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # ÂNGULOS KLEIN-NISHINA
    # --------------------------------------------------------

    @staticmethod
    def plot_angles(dataset):

        compton = dataset[
            dataset["interaction"] == "compton"
        ]

        plt.figure(figsize=(10, 6))

        plt.hist(

            compton["theta"],

            bins=80
        )

        plt.title(

            "Klein-Nishina Angular Distribution",

            fontsize=16
        )

        plt.xlabel("Angle (degrees)")
        plt.ylabel("Frequency")

        plt.savefig(

            "outputs/images/klein_nishina_angles.png",

            dpi=300,

            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # ENERGIA
    # --------------------------------------------------------

    @staticmethod
    def plot_energy_distribution(dataset):

        plt.figure(figsize=(10, 6))

        plt.hist(

            dataset["energy_after"],

            bins=80
        )

        plt.title(

            "Photon Energy Distribution",

            fontsize=16
        )

        plt.xlabel("Energy (keV)")
        plt.ylabel("Frequency")

        plt.savefig(

            "outputs/images/energy_distribution.png",

            dpi=300,

            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # TRAJETÓRIAS
    # --------------------------------------------------------

    @staticmethod
    def plot_trajectories(trajectories):

        plt.figure(figsize=(12, 12))

        for traj in trajectories[:1000]:

            plt.plot(

                traj["x"],

                traj["y"],

                alpha=0.05
            )

        plt.title(

            "Photon Trajectories",

            fontsize=18
        )

        plt.xlim(0, 400)
        plt.ylim(0, 400)

        plt.grid(alpha=0.2)

        plt.savefig(

            "outputs/images/trajectories.png",

            dpi=300,

            bbox_inches="tight"
        )

        plt.show()


# ============================================================
# ANÁLISE
# ============================================================

class Analysis:

    @staticmethod
    def generate_statistics(results):

        dataset = results["dataset"]

        stats = {

            "photoelectric_events":

                results["photoelectric"],

            "compton_events":

                results["compton"],

            "mean_energy_after":

                dataset["energy_after"].mean(),

            "mean_scattering_angle":

                dataset[
                    dataset["interaction"]
                    == "compton"
                ]["theta"].mean(),

            "total_events":

                len(dataset)
        }

        return pd.DataFrame([stats])


# ============================================================
# CONFIGURAÇÕES
# ============================================================

NUM_PARTICLES = 30000

INITIAL_ENERGY = 120.0

MATERIAL_NAME = "water"

WIDTH = 400
HEIGHT = 400

# ============================================================
# EXECUTAR
# ============================================================

material = MATERIALS[MATERIAL_NAME]

print("="*70)
print("MONTE CARLO PHOTON TRANSPORT")
print("REAL KLEIN-NISHINA MODEL")
print("="*70)

print(f"Material: {material.name}")
print(f"Initial Energy: {INITIAL_ENERGY} keV")
print(f"Particles: {NUM_PARTICLES}")

print("="*70)

# ============================================================
# SIMULAÇÃO
# ============================================================

simulation = MonteCarloSimulation(

    material=material,

    num_particles=NUM_PARTICLES,

    initial_energy=INITIAL_ENERGY,

    width=WIDTH,

    height=HEIGHT
)

results = simulation.run()

dataset = results["dataset"]

# ============================================================
# SALVAR DATASET
# ============================================================

dataset.to_csv(

    "outputs/datasets/klein_nishina_dataset.csv",

    index=False
)

# ============================================================
# ESTATÍSTICAS
# ============================================================

stats = Analysis.generate_statistics(results)

print("\nSimulation Statistics:\n")

print(stats)

# ============================================================
# VISUALIZAÇÕES
# ============================================================

Visualization.plot_dose_map(

    results["dose_map"]
)

Visualization.plot_angles(

    dataset
)

Visualization.plot_energy_distribution(

    dataset
)

Visualization.plot_trajectories(

    results["trajectories"]
)

# ============================================================
# FINAL
# ============================================================

print("\nSimulation completed successfully!")

print("\nGenerated files:")

print("✔ dose_map.png")
print("✔ klein_nishina_angles.png")
print("✔ energy_distribution.png")
print("✔ trajectories.png")
print("✔ klein_nishina_dataset.csv")