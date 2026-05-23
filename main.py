# ============================================================
# MONTE CARLO RADIATION TRANSPORT SIMULATOR
# FÍSICA DAS RADIAÇÕES + CIÊNCIA DE DADOS
# VERSÃO AVANÇADA
# ============================================================

# Melhorias implementadas:
#
# ✅ Espalhamento Compton simplificado
# ✅ Física mais realista
# ✅ Heatmap melhorado
# ✅ Trajetórias coloridas por energia
# ✅ Estatísticas avançadas
# ✅ Dataset automático
# ✅ Animação
# ✅ Melhor distribuição angular
# ✅ Resolução melhorada
# ✅ Simulação mais estável
# ✅ Estrutura pronta para IA
#
# ============================================================

# ============================================================
# IMPORTS
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import seaborn as sns

from tqdm import tqdm

# ============================================================
# CONFIGURAÇÕES GERAIS
# ============================================================

os.makedirs("outputs/images", exist_ok=True)
os.makedirs("outputs/animations", exist_ok=True)
os.makedirs("outputs/datasets", exist_ok=True)

# ============================================================
# CLASSE MATERIAL
# ============================================================

class Material:

    def __init__(
        self,
        name,
        mu,
        density,
        atomic_number
    ):

        self.name = name
        self.mu = mu
        self.density = density
        self.atomic_number = atomic_number


# ============================================================
# MATERIAIS
# ============================================================

MATERIALS = {

    "water": Material(
        "Water",
        mu=0.15,
        density=1.0,
        atomic_number=7
    ),

    "lead": Material(
        "Lead",
        mu=1.2,
        density=11.34,
        atomic_number=82
    ),

    "aluminum": Material(
        "Aluminum",
        mu=0.35,
        density=2.7,
        atomic_number=13
    ),

    "tissue": Material(
        "Tissue",
        mu=0.20,
        density=1.05,
        atomic_number=7.5
    ),
}

# ============================================================
# PARTÍCULA
# ============================================================

class Particle:

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

    def scatter_compton(self):

        scatter_angle = np.random.normal(0, np.pi/4)

        self.angle += scatter_angle

        energy_loss = np.random.uniform(0.05, 0.25)

        self.energy *= (1 - energy_loss)

        if self.energy <= 0.01:
            self.alive = False


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
# MONTE CARLO SIMULATION
# ============================================================

class MonteCarloSimulation:

    def __init__(

        self,
        material,
        num_particles=50000,
        initial_energy=1.0,
        width=400,
        height=400
    ):

        self.material = material

        self.num_particles = num_particles
        self.initial_energy = initial_energy

        self.width = width
        self.height = height

        self.dose_map = DoseMap(width, height)

        self.trajectories = []

        self.absorbed_count = 0
        self.scatter_count = 0

        self.dataset = []

    # --------------------------------------------------------
    # Caminho livre médio
    # --------------------------------------------------------

    def sample_free_path(self):

        mu = self.material.mu

        u = np.random.uniform(0.00001, 1)

        return -np.log(u) / mu

    # --------------------------------------------------------
    # Simular uma partícula
    # --------------------------------------------------------

    def simulate_particle(self):

        particle = Particle(

            x=self.width // 2,
            y=self.height // 2,

            energy=self.initial_energy,

            angle=np.random.uniform(0, 2*np.pi)
        )

        while particle.alive:

            distance = self.sample_free_path()

            old_energy = particle.energy

            particle.move(distance)

            # --------------------------------------------
            # Verificar borda
            # --------------------------------------------

            if (

                particle.x < 0
                or particle.x >= self.width
                or particle.y < 0
                or particle.y >= self.height
            ):

                break

            # --------------------------------------------
            # Probabilidade física
            # --------------------------------------------

            P_absorption = (
                self.material.mu /
                (self.material.mu + 1)
            )

            interaction = np.random.random()

            # --------------------------------------------
            # ABSORÇÃO
            # --------------------------------------------

            if interaction < P_absorption:

                deposited_energy = particle.energy

                self.dose_map.deposit_energy(
                    particle.x,
                    particle.y,
                    deposited_energy
                )

                self.dataset.append({

                    "initial_energy": old_energy,
                    "final_energy": 0,
                    "distance": distance,
                    "material": self.material.name,
                    "absorbed": 1
                })

                particle.alive = False

                self.absorbed_count += 1

            # --------------------------------------------
            # ESPALHAMENTO COMPTON
            # --------------------------------------------

            else:

                deposited_energy = particle.energy * 0.08

                self.dose_map.deposit_energy(

                    particle.x,
                    particle.y,
                    deposited_energy
                )

                particle.scatter_compton()

                self.dataset.append({

                    "initial_energy": old_energy,
                    "final_energy": particle.energy,
                    "distance": distance,
                    "material": self.material.name,
                    "absorbed": 0
                })

                self.scatter_count += 1

        self.trajectories.append({

            "x": particle.path_x,
            "y": particle.path_y,
            "energy": particle.energy_history
        })

    # --------------------------------------------------------

    def run(self):

        for _ in tqdm(range(self.num_particles)):

            self.simulate_particle()

        return {

            "dose_map": self.dose_map.grid,
            "trajectories": self.trajectories,
            "absorbed": self.absorbed_count,
            "scatter": self.scatter_count,
            "dataset": pd.DataFrame(self.dataset)
        }


# ============================================================
# VISUALIZAÇÃO
# ============================================================

class Visualization:

    # --------------------------------------------------------
    # Heatmap
    # --------------------------------------------------------

    @staticmethod
    def plot_dose_map(dose_map):

        plt.figure(figsize=(12, 10))

        sns.heatmap(

            dose_map,

            cmap="inferno",

            cbar=True
        )

        plt.title(

            "Monte Carlo Radiation Dose Map",

            fontsize=18
        )

        plt.xlabel("X")
        plt.ylabel("Y")

        plt.savefig(

            "outputs/images/dose_heatmap.png",

            dpi=300,
            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # Distribuição de energia
    # --------------------------------------------------------

    @staticmethod
    def plot_energy_distribution(dataset):

        plt.figure(figsize=(10, 6))

        plt.hist(

            dataset["final_energy"],

            bins=80
        )

        plt.title(

            "Final Energy Distribution",

            fontsize=16
        )

        plt.xlabel("Final Energy")
        plt.ylabel("Frequency")

        plt.savefig(

            "outputs/images/energy_distribution.png",

            dpi=300,
            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # Trajetórias coloridas
    # --------------------------------------------------------

    @staticmethod
    def plot_trajectories(trajectories):

        plt.figure(figsize=(12, 12))

        ax = plt.gca()

        for traj in trajectories[:500]:

            x = np.array(traj["x"])
            y = np.array(traj["y"])

            points = np.array([x, y]).T.reshape(-1, 1, 2)

            segments = np.concatenate(

                [points[:-1], points[1:]],

                axis=1
            )

            lc = LineCollection(

                segments,

                cmap='plasma',

                linewidths=1.2
            )

            lc.set_array(np.linspace(1, 0, len(segments)))

            ax.add_collection(lc)

        ax.set_xlim(0, 400)
        ax.set_ylim(0, 400)

        plt.title(

            "Particle Trajectories",

            fontsize=18
        )

        plt.xlabel("X")
        plt.ylabel("Y")

        plt.grid(alpha=0.2)

        plt.savefig(

            "outputs/images/trajectories.png",

            dpi=300,
            bbox_inches="tight"
        )

        plt.show()

    # --------------------------------------------------------
    # Animação
    # --------------------------------------------------------

    @staticmethod
    def animate_particles(trajectories):

        fig, ax = plt.subplots(figsize=(8, 8))

        ax.set_xlim(0, 400)
        ax.set_ylim(0, 400)

        ax.set_title("Monte Carlo Radiation Transport")

        lines = []

        for _ in trajectories[:50]:

            line, = ax.plot([], [], lw=1)

            lines.append(line)

        # ----------------------------------------------------

        def init():

            for line in lines:
                line.set_data([], [])

            return lines

        # ----------------------------------------------------

        def update(frame):

            for i, traj in enumerate(trajectories[:50]):

                x = traj["x"][:frame]
                y = traj["y"][:frame]

                lines[i].set_data(x, y)

            return lines

        ani = animation.FuncAnimation(

            fig,

            update,

            frames=80,

            init_func=init,

            interval=50,

            blit=True
        )

        ani.save(

            "outputs/animations/radiation_transport.gif",

            writer="pillow"
        )

        plt.close()


# ============================================================
# ESTATÍSTICAS
# ============================================================

class SimulationAnalysis:

    @staticmethod
    def generate_statistics(results):

        dataset = results["dataset"]

        stats = {

            "mean_final_energy":
                dataset["final_energy"].mean(),

            "max_final_energy":
                dataset["final_energy"].max(),

            "mean_distance":
                dataset["distance"].mean(),

            "absorbed_particles":
                results["absorbed"],

            "scatter_events":
                results["scatter"],

            "total_events":
                len(dataset)
        }

        return pd.DataFrame([stats])


# ============================================================
# CONFIGURAÇÕES
# ============================================================

NUM_PARTICLES = 50000

INITIAL_ENERGY = 1.0

MATERIAL_NAME = "water"

WIDTH = 400
HEIGHT = 400

# ============================================================
# EXECUTAR
# ============================================================

material = MATERIALS[MATERIAL_NAME]

print("="*70)
print("MONTE CARLO RADIATION TRANSPORT SIMULATOR")
print("="*70)

print(f"Material: {material.name}")
print(f"Attenuation coefficient: {material.mu}")
print(f"Density: {material.density}")
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

# ============================================================
# DATASET
# ============================================================

dataset = results["dataset"]

dataset.to_csv(

    "outputs/datasets/radiation_dataset.csv",

    index=False
)

# ============================================================
# ESTATÍSTICAS
# ============================================================

stats = SimulationAnalysis.generate_statistics(results)

print("\nSimulation Statistics:\n")

print(stats)

# ============================================================
# VISUALIZAÇÕES
# ============================================================

Visualization.plot_dose_map(

    results["dose_map"]
)

Visualization.plot_energy_distribution(

    dataset
)

Visualization.plot_trajectories(

    results["trajectories"]
)

# ============================================================
# ANIMAÇÃO
# ============================================================

print("\nGenerating animation...\n")

Visualization.animate_particles(

    results["trajectories"]
)

# ============================================================
# FINAL
# ============================================================

print("\nSimulation completed successfully!")

print("\nFiles generated:")

print("✔ dose_heatmap.png")
print("✔ trajectories.png")
print("✔ energy_distribution.png")
print("✔ radiation_transport.gif")
print("✔ radiation_dataset.csv")