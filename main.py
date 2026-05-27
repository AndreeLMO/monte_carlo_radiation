# ============================================================
# MONTE CARLO RADIATION TRANSPORT
# KLEIN-NISHINA REAL
# ============================================================
# Autor: André Luiz Magalhães de Oliveira
# Formação: Físico Médico | Especialista em Data Science & Analytics
# Universidade de São Paulo
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# ============================================================
# CONFIGURAÇÕES
# ============================================================

os.makedirs("outputs/images", exist_ok=True)
os.makedirs("outputs/datasets", exist_ok=True)

ELECTRON_REST_ENERGY = 511.0  # keV
CLASSICAL_ELECTRON_RADIUS = 2.817e-13

# ============================================================
# MATERIAL
# ============================================================

class Material:
    def __init__(self, name, density, atomic_number, mu_total):
        self.name = name
        self.density = density
        self.atomic_number = atomic_number
        self.mu_total = mu_total

MATERIALS = {
    "water": Material("Water", 1.0, 7, 0.15),
    "tissue": Material("Tissue", 1.05, 7.4, 0.20),
    "aluminum": Material("Aluminum", 2.7, 13, 0.35),
    "lead": Material("Lead", 11.34, 82, 1.20)
}

# ============================================================
# PROBABILIDADES FÍSICAS
# ============================================================

def calculate_interaction_probabilities(energy, material):
    Z = material.atomic_number
    photoelectric = Z**4 / (energy + 1e-6)**3
    compton = Z / (energy + 1e-6)
    total = photoelectric + compton
    return photoelectric/total, compton/total

def klein_nishina_differential_cs(energy, theta):
    epsilon = energy / ELECTRON_REST_ENERGY
    E_ratio = 1 / (1 + epsilon * (1 - np.cos(theta)))
    term1 = E_ratio**2
    term2 = E_ratio + 1/E_ratio - np.sin(theta)**2
    return (CLASSICAL_ELECTRON_RADIUS**2 / 2) * term1 * term2

def sample_klein_nishina_angle(energy, n_samples=1000):
    angles = np.linspace(0, np.pi, n_samples)
    probs = np.array([klein_nishina_differential_cs(energy, t) for t in angles])
    probs /= probs.sum()
    return np.random.choice(angles, p=probs)

# ============================================================
# FÓTON
# ============================================================

class Photon:
    def __init__(self, x, y, energy, angle):
        self.x, self.y, self.energy, self.angle = x, y, energy, angle
        self.alive = True
        self.path_x, self.path_y = [x], [y]
        self.energy_history = [energy]

    def move(self, distance):
        self.x += distance * np.cos(self.angle)
        self.y += distance * np.sin(self.angle)
        self.path_x.append(self.x); self.path_y.append(self.y)
        self.energy_history.append(self.energy)

    def compton_scatter(self):
        theta = sample_klein_nishina_angle(self.energy)
        E_initial = self.energy
        E_scattered = E_initial / (1 + (E_initial/ELECTRON_REST_ENERGY) * (1 - np.cos(theta)))
        deposited = E_initial - E_scattered
        self.energy, self.angle = E_scattered, self.angle + theta
        if self.energy < 1: self.alive = False
        return deposited, theta

    def photoelectric_absorption(self):
        deposited = self.energy
        self.energy, self.alive = 0, False
        return deposited

# ============================================================
# DOSE MAP
# ============================================================

class DoseMap:
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.grid = np.zeros((height, width))

    def deposit_energy(self, x, y, energy):
        ix, iy = int(x), int(y)
        if 0 <= ix < self.width and 0 <= iy < self.height:
            self.grid[iy, ix] += energy

# ============================================================
# SIMULAÇÃO MONTE CARLO
# ============================================================

class MonteCarloSimulation:
    def __init__(self, material, num_particles=30000, initial_energy=120.0, width=400, height=400):
        self.material, self.num_particles, self.initial_energy = material, num_particles, initial_energy
        self.width, self.height = width, height
        self.dose_map, self.dataset, self.trajectories = DoseMap(width, height), [], []
        self.photoelectric_events, self.compton_events = 0, 0

    def sample_free_path(self):
        return -np.log(np.random.uniform(1e-10, 1)) / self.material.mu_total

    def simulate_photon(self):
        photon = Photon(self.width//2, self.height//2, self.initial_energy, np.random.uniform(0, 2*np.pi))
        while photon.alive:
            distance = self.sample_free_path()
            photon.move(distance)
            if not (0 <= photon.x < self.width and 0 <= photon.y < self.height): break
            P_photo, P_comp = calculate_interaction_probabilities(photon.energy, self.material)
            if np.random.random() < P_photo:
                deposited = photon.photoelectric_absorption()
                self.dose_map.deposit_energy(photon.x, photon.y, deposited)
                self.photoelectric_events += 1
                self.dataset.append({"interaction":"photoelectric","energy_before":deposited,"energy_after":0,"theta":0,"distance":distance,"x":photon.x,"y":photon.y})
            else:
                energy_before = photon.energy
                deposited, theta = photon.compton_scatter()
                self.dose_map.deposit_energy(photon.x, photon.y, deposited)
                self.compton_events += 1
                self.dataset.append({"interaction":"compton","energy_before":energy_before,"energy_after":photon.energy,"theta":np.degrees(theta),"distance":distance,"x":photon.x,"y":photon.y})
        self.trajectories.append({"x":photon.path_x,"y":photon.path_y})

    def run(self):
        for _ in tqdm(range(self.num_particles)): self.simulate_photon()
        return {"dose_map":self.dose_map.grid,"dataset":pd.DataFrame(self.dataset),"trajectories":self.trajectories,"photoelectric":self.photoelectric_events,"compton":self.compton_events}

# ============================================================
# ANÁLISE
# ============================================================

class Analysis:
    @staticmethod
    def generate_statistics(results):
        dataset = results["dataset"]
        return pd.DataFrame([{
            "photoelectric_events": results["photoelectric"],
            "compton_events": results["compton"],
            "mean_energy_after": dataset["energy_after"].mean(),
            "mean_scattering_angle": dataset[dataset["interaction"]=="compton"]["theta"].mean(),
            "total_events": len(dataset)
        }])

# ============================================================
# EXECUÇÃO
# ============================================================

NUM_PARTICLES, INITIAL_ENERGY, MATERIAL_NAME = 30000, 120.0, "water"
material = MATERIALS[MATERIAL_NAME]

print("="*70)
print("MONTE CARLO PHOTON TRANSPORT — REAL KLEIN-NISHINA MODEL")
print("="*70)
print(f"Material: {material.name} | Initial Energy: {INITIAL_ENERGY} keV | Particles: {NUM_PARTICLES}")
print("="*70)

simulation = MonteCarloSimulation(material, NUM_PARTICLES, INITIAL_ENERGY)
results = simulation.run()
dataset = results["dataset"]

dataset.to_csv("outputs/datasets/klein_nishina_dataset.csv", index=False)
print("\nSimulation Statistics:\n", Analysis.generate_statistics(results))
print("\nSimulation completed successfully!")
