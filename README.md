# ☢️ Simulação de Transporte de Radiação Ionizante via Método de Monte Carlo

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Physics](https://img.shields.io/badge/Physics-Medical%20%26%20Nuclear-orange.svg)]()
[![Framework](https://img.shields.io/badge/Framework-OOP%20Architecture-green.svg)]()

Este repositório aloja um simulador computacional estocástico bidimensional e tridimensional, desenvolvido em Python, concebido para modelar o transporte, espalhamento e atenuação de fótons de alta energia (radiação eletromagnética ionizante) ao atravessarem meios materiais homogéneos. 

Utilizando o **Método de Monte Carlo (MMC)**, o sistema reconstrói o histórico individual (histórias) de até $10^6$ partículas, amostrando variáveis probabilísticas macroscópicas a partir de seções de choque microestruturais. O objetivo principal é mapear com precisão clínica a deposição espacial de energia (**Dose Absorvida**) em estruturas digitais (*phantoms*).

---

## 📚 1. Fundamentação Teórica e Física das Radiações

O transporte de radiação ionizante na matéria é governado por processos estocásticos de colisão e transferência de energia cinética. Matematicamente, a atenuação de um feixe de fótons monoenergéticos colimados é modelada pela equação diferencial de transporte unidimensional:

$$\frac{dN}{dx} = -\mu \cdot N(x)$$

A integração desta equação ao longo de uma espessura finita $x$ resulta na **Lei de Beer-Lambert**:

$$N(x) = N_0 \cdot e^{-\mu \cdot x}$$

Onde:
* $N_0$ representa o número de fótons incidentes no ponto de entrada do meio.
* $x$ é a espessura linear penetrada no material (expressa em $\text{cm}$).
* $\mu$ é o **Coeficiente de Atenuação Linear Total** ($\text{cm}^{-1}$), que define a probabilidade de interação por unidade de comprimento.

### 🎲 Amostragem do Livre Caminho Médio (LCM)
O percurso ou distância livre $s$ que um fóton percorre entre duas interações consecutivas é uma variável aleatória cuja Função Densidade de Probabilidade (FDP) é dada por:

$$f(s) = \mu \cdot e^{-\mu \cdot s} \quad \text{para} \quad s \ge 0$$

Para simular computacionalmente esta distância, aplica-se o **Método da Transformação Inversa**. Igualando a Função Distribuição Cumulativa (FDC) a um número pseudoaleatório $u$, uniformemente distribuído no intervalo estatístico $(0, 1]$, obtém-se:

$$u = F(s) = \int_{0}^{s} \mu \cdot e^{-\mu \cdot t} \, dt = 1 - e^{-\mu \cdot s}$$

Isolando-se matematicamente a variável de interesse $s$, e considerando que a distribuição de $1 - u$ é estatisticamente idêntica à de $u$, chegamos à equação operacional do simulador:

$$s = -\frac{\ln(u)}{\mu}$$

### ⚡ Fenómenos de Interação Modelados
* **Absorção Fotoelétrica ($\tau$):** O fóton colide com um elétron orbital fortemente ligado, transferindo-lhe $100\%$ da sua energia cinética. O fóton é extinto do sistema (fim da história da partícula), e a energia é depositada integralmente no pixel da colisão.
* **Espalhamento Compton ($\sigma$):** O fóton interage com um elétron fracamente ligado ou livre. Ocorre uma colisão elástica onde o fóton é defletido por um ângulo $\theta$ em relação à sua direção original, sofrendo uma perda parcial de energia proporcional à severidade do desvio angular, modelada no algoritmo por uma taxa de atenuação fracionária estável de $\Delta E = 15\%$.

---

## 🛠️ 2. Arquitetura do Software e Engenharia de Código

O simulador foi estruturado seguindo o paradigma de **Programação Orientada a Objetos (POO)**, maximizando a modularidade, escalabilidade e a analogia física direta no código.
