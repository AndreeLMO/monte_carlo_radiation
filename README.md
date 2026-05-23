# Simulação de Monte Carlo para Transporte de Radiação: Modelo Klein-Nishina Real

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Physics](https://img.shields.io/badge/physics-Medical%20%7C%20Radiological-red.svg)]()

Este repositório hospeda o desenvolvimento de um simulador estocástico nativo em Python projetado para modelar o transporte, espalhamento e deposição de energia (dose) de fótons na faixa de energias de radiodiagnóstico e tomografia computadorizada ($10 \text{ keV}$ a $150 \text{ keV}$) em meios condensados homogêneos. O núcleo analítico do algoritmo realiza a amostragem estocástica exata da **Seção de Choque Diferencial de Klein-Nishina** através da inversão numérica da Função de Distribuição Acumulada (FDA), superando aproximações analíticas simplificadas (como Henyey-Greenstein) comumente aplicadas na literatura.

---

## 📋 Sumário
1. [🔬 Fundamentação Teórica Expandida](#-fundamentação-teórica-expandida)
    * [Cinemática do Transporte Macro-Micro](#cinemática-do-transporte-macro-micro)
    * [O Efeito Fotoelétrico e Limites de Energia](#o-efeito-fotoelétrico-e-limites-de-energia)
    * [Dedução Mecânico-Quântica do Espalhamento Compton](#dedução-mecânico-quântica-do-espalhamento-compton)
    * [A Formulação de Klein-Nishina](#a-formulação-de-klein-nishina)
2. [💻 Implementação Algorítmica](#-implementação-algorítmica)
    * [Amostragem de Caminhos e Geometria](#amostragem-de-caminhos-e-geometria)
3. [📊 Análise Avançada de Resultados e Imagens](#-análise-avançada-de-resultados-e-imagens)
    * [Consolidação Numérica de Saída](#consolidação-numérica-de-saída)
    * [Figura 1: Distribuição Espacial de Dose (Isodose)](#figura-1-distribuição-espacial-de-dose-isodose)
    * [Figura 2: Validação Estatística Angular](#figura-2-validação-estatística-angular)
    * [Figura 3: Espectroscopia de Fótons e Contínuo Compton](#figura-3-espectroscopia-de-fótons-e-contínuo-compton)
4. [🛠️ Guia de Execução](#%EF%B8%8F-guia-de-execução)

---

## 🔬 Fundamentação Teórica Expandida

### Cinemática do Transporte Macro-Micro
O deslocamento de um fóton em um meio material condensado é um processo puramente probabilístico descrito pela Teoria do Transporte de Partículas. O coeficiente de atenuação linear total $\mu(E, Z)$ define a probabilidade de uma partícula sofrer *qualquer* tipo de colisão por unidade de comprimento (geralmente expresso em $\text{cm}^{-1}$). 

Quando um feixe com $N_0$ fótons incide em uma espessura $x$, o número de fótons que atravessam sem interagir segue a Lei de Beer-Lambert:
$$N(x) = N_0 \cdot e^{-\mu x}$$

No nível microscópico de simulação individual (Histórico de Monte Carlo), a distância linear $s$ percorrida por um fóton entre dois vértices consecutivos de interação (Livre Caminho Médio Estocástico) é calculada igualando a probabilidade acumulada a um número pseudoaleatório $U$ uniformemente distribuído no intervalo $(0,1)$:
$$s = -\frac{\ln(U)}{\mu(E, Z)}$$

### O Efeito Fotoelétrico e Limites de Energia
O Efeito Fotoelétrico caracteriza a absorção total do fóton por um elétron ligado a uma das camadas orbitais internas do átomo (predominantemente a camada K). A energia do fóton incidente $E$ é completamente transferida: uma fração supera a energia de ligação do elétron ($B_e$) e o restante é convertido em energia cinética do fotoelétron ejetado ($E_c = E - B_e$).

A seção de choque atômica para este processo ($\tau$) varia drasticamente com o número atômico do absorvedor ($Z$) e com a energia cinética do fóton:
$$\tau \approx \text{Constante} \cdot \frac{Z^4}{E^3}$$



Devido à dependência com $Z^4$, o efeito fotoelétrico é o mecanismo primário de contraste em radiologia médica (diferenciando estruturas ósseas de tecidos moles) e eficiência de blindagens radiológicas. Contudo, em meios de baixo número atômico como a água ($Z_{\text{ef}} \approx 7.4$), a probabilidade de absorção fotoelétrica decai exponencialmente para energias superiores a $50 \text{ keV}$, cedendo a dominância estatística ao efeito Compton.

### Dedução Mecânico-Quântica do Espalhamento Compton
O espalhamento Compton ocorre quando um fóton colide com um elétron considerado livre e em repouso (elétrons de valência cujas energias de ligação são desprezíveis frente à energia do fóton incidente). 



Ao aplicar a conservação do momentum linear relativístico $\vec{p}$ e da energia total $E$ no sistema bidimensional formado pelo fóton defletido em um ângulo polar $\theta$ e pelo elétron de recuo projetado em um ângulo $\psi$, obtém-se a consagrada relação de degradação do comprimento de onda de Compton ($\lambda' - \lambda$):
$$\lambda' - \lambda = \frac{h}{m_e c}(1 - \cos\theta)$$

Convertendo comprimentos de onda para unidades de energia ($E = hc/\lambda$), a energia do fóton pós-espalhamento ($E'$) assume a forma matemática:
$$E' = \frac{E}{1 + \epsilon(1 - \cos\theta)}$$

Onde $\epsilon = \frac{E}{m_e c^2}$ representa a energia reduzida do fóton em unidades de energia de repouso do elétron ($m_e c^2 = 511.0 \text{ keV}$). A energia cinética transferida ao elétron de recuo, que constitui a dose depositada localmente no voxel da colisão, é dada por:
$$T_e = E - E' = E \left[ \frac{\epsilon(1 - \cos\theta)}{1 + \epsilon(1 - \cos\theta)} \right]$$

### A Formulation de Klein-Nishina
Embora a cinemática angular seja determinística pelas leis de conservação, a probabilidade física de o fóton ser defletido em um ângulo específico $\theta$ por unidade de ângulo sólido ($d\Omega$) exige o tratamento quântico relativístico da equação de Dirac. Desenvolvida por Oskar Klein e Yoshio Nishina, a seção de choque diferencial por elétron livre é expressa como:

$$\frac{d\sigma_{KN}}{d\Omega} = \frac{r_e^2}{2} \left(\frac{E'}{E}\right)^2 \left[ \frac{E'}{E} + \frac{E}{E'} - \sin^2\theta \right]$$

Onde $r_e = 2.817 \times 10^{-13} \text{ cm}$ é o raio clássico do elétron. À medida que a energia incidente $E$ cresce, a distribuição angular perde sua característica de simetria simétrica (fórmula de espalhamento Thomson clássica) e projeta-se fortemente em direção frontal (ângulos agudos).

---

## 💻 Implementação Algorítmica

O fluxograma operacional do código segue o rastreamento individual de históricos até que critérios de corte geométricos ou energéticos sejam satisfeitos.

## 💻 Fluxograma Operacional (Monte Carlo)

```mermaid
graph TD
    A([Início: Injeção de Fóton Primário<br>E0, x0, y0]) --> B[Calcular Coeficientes Totais &mu;]
    B --> C[Sorteio do Passo<br>s = -ln U / &mu;]
    C --> D[Atualizar Posição Espacial]
    D --> E{Fora do Alvo?}
    
    E -- Sim --> F([Fim do Histórico:<br>Próxima Partícula])
    E -- Não --> G[Sorteio do Tipo de Interação]
    
    G --> H{Qual Interação?}
    H -- Fotoelétrico --> I[Absorção Total]
    I --> J([Morte da Partícula])
    J --> F
    
    H -- Compton --> K[Amostragem Angular de Klein-Nishina &theta;<br>Sorteio Azimutal Uniforme &phi; = 2&pi;U]
    K --> L[Calcular Nova Energia E' e Tr<br>Acumular Deposição de Dose Local Tr]
    L --> M{Critério de Corte:<br>E' < 1 keV?}
    
    M -- Sim --> J
    M -- Não --> N[Loop: Atualizar E = E']
    N --> B
    
    %% Estilização para o gráfico ficar elegante no GitHub
    style A fill:#64b5f6,stroke:#1565c0,stroke-width:2px,color:#000
    style F fill:#90caf9,stroke:#1565c0,stroke-width:2px,color:#000
    style J fill:#ef9a9a,stroke:#c62828,stroke-width:2px,color:#000
    style E fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style H fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style M fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
---

## 4. Análise Avançada de Resultados e Discussão

A simulação de validação padrão adotou um modelo estocástico rigoroso lançando **$30.000$ fótons primários** monoenergéticos com energia inicial de **$120,0 \text{ keV}$**. O ponto de injeção do feixe foi configurado estritamente no centro geométrico da matriz de transporte, correspondendo às coordenadas $(x_c, y_c) = (200 \text{ mm}, 200 \text{ mm})$, inserido num simulador cúbico (phantom) homogêneo de $400 \times 400 \text{ mm}^2$. 

Abaixo, detalha-se a dinâmica operacional do algoritmo e a resposta físico-estatística do meio (Água líquida pura, $Z_{\text{ef}} \approx 7.4$).

---

### 4.1. Fluxograma Operacional do Histórico de Partículas

Para garantir a estabilidade e a reprodução fiel da física local sem perda de alinhamento textual no GitHub, o mapeamento lógico do transporte de cada histórico foi estruturado através da sintaxe nativa `mermaid`:

```mermaid
graph TD
    A([Início: Injeção de Fóton Primário<br>E0 = 120 keV, x0=200, y0=200]) --> B[Calcular Coeficientes Totais &mu; E, Z]
    B --> C[Sorteio do Passo Estocástico<br>s = -ln U / &mu;]
    C --> D[Atualizar Posição Cartesiana]
    D --> E{Fora dos Limites<br>do Alvo?}
    
    E -- Sim --> F([Fim do Histórico:<br>Contabilizar Fuga])
    E -- Não --> G[Sorteio do Canal de Interação]
    
    G --> H{Qual o Evento<br>Sorteado?}
    H -- Fotoelétrico --> I[Absorção Total da Energia]
    I --> J([Morte da Partícula])
    J --> F
    
    H -- Compton --> K[Amostragem Angular Real de Klein-Nishina &theta;<br>Sorteio Azimutal Uniforme &phi; = 2&pi;U]
    K --> L[Calcular Nova Energia E' e Tr<br>Acumular Deposição de Dose Local Tr]
    L --> M{Critério de Corte<br>Energético: E' < 1 keV?}
    
    M -- Sim --> J
    M -- Não --> N[Loop: Atualizar E = E']
    N --> B
    
    style A fill:#64b5f6,stroke:#1565c0,stroke-width:2px,color:#000
    style F fill:#90caf9,stroke:#1565c0,stroke-width:2px,color:#000
    style J fill:#ef9a9a,stroke:#c62828,stroke-width:2px,color:#000
    style E fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style H fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
    style M fill:#fff59d,stroke:#fbc02d,stroke-width:2px,color:#000
