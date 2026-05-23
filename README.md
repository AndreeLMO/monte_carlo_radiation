# Simulação de Monte Carlo para Transporte de Radiação: Modelo Klein-Nishina Real

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Contributions](https://img.shields.io/badge/contributions-welcome-orange.svg)](CONTRIBUTING.md)

Este repositório contém o código-fonte e a documentação técnica de um simulador estocástico baseado no **Método de Monte Carlo** para o transporte e deposição de dose de fótons de baixa e média energia (faixa de radiodiagnóstico e tomografia computadorizada) em meios condensados homogêneos. O diferencial deste código é a amostragem analítica exata da **Seção de Choque Diferencial de Klein-Nishina**, sem aproximações numéricas grosseiras.

---

## 📋 Sumário
1. [Fundamentação Teórica](#-fundamentação-teórica)
    * [Coeficiente de Atenuação Linear](#coeficiente-de-atenuação-linear)
    * [Efeito Fotoelétrico](#efeito-fotoelétrico)
    * [Espalhamento Compton Incoerente](#espalhamento-compton-incoerente)
    * [Seção de Choque de Klein-Nishina](#seção-de-choque-de-klein-nishina)
2. [Arquitetura do Código](#-arquitetura-do-código)
    * [Modelagem de Materiais](#modelagem-de-materiais)
    * [Amostragem Estocástica Angular](#amostragem-estocástica-angular)
3. [Resultados e Análise Físico-Estatística](#-resultados-e-análise-físico-estatística)
    * [Estatística de Interações](#estatística-de-interações)
    * [Perfil Espacial de Deposição de Dose](#perfil-espacial-de-deposição-de-dose)
    * [Distribuição Angular e Degradação Energética](#distribuição-angular-e-degradação-energética)
4. [Como Executar o Projeto](#-como-executar-o-projeto)

---

## 🔬 Fundamentação Teórica

### Coeficiente de Atenuação Linear
O transporte macroscópico de fótons através da matéria é probabilisticamente governado pelo coeficiente de atenuação linear total $\mu(E, Z)$, que representa a probabilidade de interação por unidade de comprimento. Na faixa energética avaliada ($10 \text{ keV}$ a $150 \text{ keV}$), as interações dominantes são expressas por:

$$\mu(E, Z) = \tau(E, Z) + \sigma(E, Z)$$

Onde $\tau$ representa o coeficiente para o efeito fotoelétrico e $\sigma$ denota o coeficiente para o espalhamento Compton.

### Efeito Fotoelétrico
O efeito fotoelétrico é um processo de absorção total, onde o fóton incidente transfere toda a sua energia cinética para um elétron fortemente ligado das camadas mais internas do átomo (camadas $K$ ou $L$). A probabilidade de ocorrência por átomo é altamente dependente do número atômico ($Z$) do meio e inversamente proporcional à energia ($E$):

$$\tau \propto \frac{Z^4}{E^3}$$

### Espalhamento Compton Incoerente
No espalhamento Compton, o fóton colide inelasticamente com um elétron periférico (considerado livre e em repouso). Aplicando as leis de conservação de momento linear quádruplo, a energia do fóton defletido $E'$ é dada pela clássica equação de Compton:

$$E' = \frac{E}{1 + \frac{E}{m_e c^2}(1 - \cos\theta)}$$

Onde $m_e c^2 = 511.0 \text{ keV}$ é a energia de repouso do elétron e $\theta$ é o ângulo polar de espalhamento.

### Seção de Choque de Klein-Nishina
A distribuição probabilística do ângulo polar $\theta$ é rigorosamente descrita pela **Seção de Choque Diferencial de Klein-Nishina (SDKN)**, derivada a partir da eletrodinâmica quântica relativística utilizando a equação de Dirac:

$$\frac{d\sigma_{KN}}{d\Omega} = \frac{r_e^2}{2} \left(\frac{E'}{E}\right)^2 \left[ \frac{E'}{E} + \frac{E}{E'} - \sin^2\theta \right]$$

Onde $r_e = 2.817 \times 10^{-13} \text{ cm}$ é o raio clássico do elétron. O ângulo azimutal $\phi$ possui simetria cilíndrica e é distribuído uniformemente entre $0$ e $2\pi$.

---

## 💻 Arquitetura do Código

O simulador foi construído em Python utilizando programação orientada a objetos (POO) combinada com processamento matemático vetorizado via `numpy`.

### Modelagem de Materiais
A classe `Material` encapsula as propriedades físicas dos alvos simulados obtidos a partir das bases de dados do NIST:

```python
class Material:
    def __init__(self, name, density, atomic_number, mu_total):
        self.name = name
        self.density = density                # g/cm³
        self.atomic_number = atomic_number    # Z efetivo
        self.mu_total = mu_total              # cm⁻¹
