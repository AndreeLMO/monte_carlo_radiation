# Simulação de Monte Carlo para Transporte de Radiação: Modelo Klein-Nishina Real

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![Physics](https://img.shields.io/badge/physics-Medical%20%7C%20Radiological-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 1. Introdução

O transporte de partículas ionizantes em meios condensados constitui o pilar fundamental da física médica, radioproteção e engenharia nuclear. Na faixa de energias correspondente ao radiodiagnóstico e à tomografia computadorizada ($10 \text{ keV}$ a $150 \text{ keV}$), a modelagem matemática exata da deposição de dose e da radiação espalhada secundária é crucial para a otimização de blindagens e cálculos dosimétricos clínicos.

Os métodos analíticos tradicionais para resolver a Equação de Transporte de Boltzmann apresentam limitações severas ao lidar com geometrias complexas e meios heterogêneos. Como alternativa, o **Método de Monte Carlo** consolidou-se como o padrão-ouro (*gold standard*) na física das radiações. Trata-se de uma abordagem estocástica onde o histórico de milhares de partículas individuais é rastreado probabilisticamente através da amostragem de funções de distribuição de probabilidade baseadas em seções de choque microscópicas reais.

Nesta faixa de energia diagnóstica, dois efeitos competitivos dominam a atenuação de fótons: o **Efeito Fotoelétrico** e o **Espalhamento Compton**. Enquanto muitas simulações simplificadas utilizam aproximações analíticas para a distribuição angular (como a aproximação de Henyey-Greenstein), este projeto implementa um simulador nativo em Python que realiza a amostragem exata da **Seção de Choque Diferencial de Klein-Nishina**, capturando a anisotropia real do espalhamento quântico-relativístico.

---

## 2. Fundamentação Teórica

### 2.1. Cinemática do Transporte Macro-Micro
A probabilidade de um fóton interagir com o meio material por unidade de comprimento é quantificada pelo coeficiente de atenuação linear total, $\mu_t(E, Z)$, expresso em $\text{cm}^{-1}$. Matematicamente, a atenuação macroscópica de um feixe monoenergético com $N_0$ fótons ao atravessar uma espessura espessa $x$ segue a Lei de Beer-Lambert:

$$N(x) = N_0 \cdot e^{-\mu_t(E, Z) x}$$

Em uma simulação de histórico individual (Monte Carlo), o livre caminho caminho percorrido por um fóton entre dois vértices consecutivos de colisão é uma variável aleatória contínua. Igualando a função de distribuição acumulada a um número pseudoaleatório $U$ uniformemente distribuído no intervalo $(0, 1]$, deduz-se a distância linear do passo estocástico ($s$):

$$s = -\frac{\ln(U)}{\mu_t(E, Z)}$$

### 2.2. O Efeito Fotoelétrico e a Absorção Total
O Efeito Fotoelétrico descreve o processo no qual um fóton incidente colide com um elétron fortemente ligado às camadas internas do átomo (predominantemente a camada K). O fóton transfere integralmente sua energia para o elétron, cessando sua existência. O fotoelétron é ejetado com energia cinética dada por:

$$E_c = E - B_e$$

Onde $B_e$ é a energia de ligação orbital. A seção de choque atômica fotoelétrica ($\tau$) possui uma dependência crítica com o número atômico ($Z$) do meio absorvedor e com a energia cinética do fóton incidente ($E$):

$$\tau \propto \frac{Z^4}{E^3}$$

Devido a essa dependência de quarta potência com $Z$, o efeito fotoelétrico é o mecanismo primário responsável pelo contraste em imagens radiográficas (distinguindo estruturas ósseas de tecidos moles) e pela alta eficiência de materiais densos como o chumbo em blindagens físicas.

### 2.3. O Espalhamento Compton e a Formulação de Klein-Nishina
O espalhamento Compton ocorre quando o fóton interage com um elétron considerado livre e em repouso (elétrons de valência com energia de ligação desprezível frente à energia do fóton). Ao aplicar a conservação do quadrimomentum relativístico no sistema bidimensional formado pela deflexão do fóton em um ângulo polar $\theta$ e o elétron de recuo em um ângulo $\psi$, obtém-se a relação cinemática para a energia do fóton pós-colisão ($E'$):

$$E' = \frac{E}{1 + \epsilon(1 - \cos\theta)}$$

Onde $\epsilon = \frac{E}{m_0 c^2}$ representa a energia reduzida do fóton em unidades da energia de repouso do elétron ($m_0 c^2 = 511,0 \text{ keV}$). A energia mecânica restante é transferida ao elétron sob a forma de energia cinética de recuo ($T_e$), que constitui a dose depositada localmente no meio absorvedor:

$$T_e = E - E' = E \left[ \frac{\epsilon(1 - \cos\theta)}{1 + \epsilon(1 - \cos\theta)} \right]$$

A probabilidade física de o fóton ser defletido em um ângulo polar $\theta$ por unidade de ângulo sólido ($d\Omega$) é governada pela **Seção de Choque Diferencial de Klein-Nishina**, derivada através da mecânica quântica relativística da equação de Dirac:

$$\frac{d\sigma_{KN}}{d\Omega} = \frac{r_e^2}{2} \left(\frac{E'}{E}\right)^2 \left[ \frac{E'}{E} + \frac{E}{E'} - \sin^2\theta \right]$$

Onde $r_e = 2,817 \times 10^{-13} \text{ cm}$ é o raio clássico do elétron. À medida que a energia incidente $E$ se eleva, a distribuição angular perde sua característica de simetria simétrica (regida pelo espalhamento Thomson clássico) e projeta-se predominantemente em direção frontal (ângulos agudos).

---

## 3. Metodologia Computacional

### 3.1. Arquitetura do Simulador e Design de Software
O simulador foi integralmente desenvolvido em Python 3.8+ utilizando o paradigma de Programação Orientada a Objetos (POO) combinada com computação vetorizada (`NumPy`) para otimização do processamento estocástico. A arquitetura divide-se nas seguintes entidades estruturais:

* **`Material`**: Classe responsável por encapsular as propriedades físico-químicas do meio atenuador, incluindo densidade ($\rho$), número atômico ($Z$) e o coeficiente de atenuação total ($\mu_t$). O simulador conta com um banco de dados calibrado para quatro meios de interesse: Água ($H_2O$), Tecido Humano, Alumínio ($Al$) e Chumbo ($Pb$).
* **`Photon`**: Classe que rastreia dinamicamente o estado físico de cada partícula individual. Controla os atributos de posição cartesiana bidimensional ($x, y$), energia instantânea ($E$), direção angular azimutal ($\phi$) e o histórico completo de seu vetor de trajetória e degradação energética.
* **`DoseMap`**: Matriz bidimensional discreta acoplada que atua como um fantoma digital numérico, registrando e acumulando espacialmente a energia cinética cedida pelo fóton ($T_e$) a cada colisão nas coordenadas espaciais correspondentes.
* **`MonteCarloSimulation`**: Motor controlador do loop estocástico. Gerencia o lançamento dos históricos de partículas, executa os sorteios probabilísticos e consolida a saída de dados brutos e estatísticos.

### 3.2. Amostragem Numérica da Distribuição de Klein-Nishina
Como a Função de Distribuição Acumulada (FDA) obtida a partir da integração da equação de Klein-Nishina não possui uma forma analítica inversível diretamente por métodos simples, o algoritmo resolve a amostragem angular de forma exata via **Inversão Numérica Discreta**. 

A seção de choque diferencial é discretizada em um espaço de alta resolução ($1.000$ nós no intervalo $[0, \pi]$). A distribuição de probabilidade cumulativa é normalizada e mapeada numericamente, permitindo que o ângulo polar $\theta$ seja sorteado stocasticamente com base no perfil quântico real de Klein-Nishina para a energia exata daquela colisão. O ângulo azimutal ($\phi$) é tratado de forma isotrópica e amostrado uniformemente:

$$\phi = 2\pi \cdot U_2, \quad U_2 \sim \mathcal{U}(0, 1)$$

### 3.3. Algoritmo de Transporte e Critérios de Parada
O fluxo operacional do transporte de cada fóton individual segue estritamente a cadeia de decisões lógicas mapeada abaixo:

```mermaid
graph TD
    A([Início: Injeção de Fóton Primário<br>E0, x0, y0]) --> B[Calcular Coeficientes Totais &mu; E, Z]
    B --> C[Sorteio do Passo Estocástico<br>s = -ln U / &mu;]
    C --> D[Atualizar Posição Espacial Cartesiana]
    D --> E{Fora dos Limites<br>do Fantoma?}
    
    E -- Sim --> F([Fim do Histórico:<br>Contabilizar Fuga])
    E -- Não --> G[Sorteio Probabilístico do Canal de Interação]
    
    G --> H{Qual Interação<br>Sorteada?}
    H -- Fotoelétrico --> I[Absorção Total da Energia]
    I --> J([Morte da Partícula])
    J --> F
    
    H -- Compton --> K[Amostragem Angular Numérica de Klein-Nishina &theta;<br>Sorteio Azimutal Uniforme &phi; = 2&pi;U]
    K --> L[Calcular Nova Energia E' e Tr<br>Acumular Deposição de Dose Local no DoseMap]
    L --> M{Critério de Corte Energético:<br>E' < 1 keV?}
    
    M -- Sim --> J
    M -- Não --> N[Loop: Atualizar E = E' para Próxima Colisão]
    N --> B
    
    style A fill:#64b5f6,stroke:#1565c0,stroke-width:2px,color:#000
    style F fill:#90caf9,stroke:#1565c0,stroke-width:2px,color:#000
    style J fill:#ef9a9a,stroke:#c62828,stroke-width:2px,color:#000
    style E fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style H fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style M fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
