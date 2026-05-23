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
3. [📁 Estrutura de Pastas e Inserção de Imagens](#-estrutura-de-pastas-e-inserção-de-imagens)
4. [📊 Análise Avançada de Resultados e Imagens](#-análise-avançada-de-resultados-e-imagens)
    * [Consolidação Numérica de Saída](#consolidação-numérica-de-saída)
    * [Figura 1: Distribuição Espacial de Dose (Isodose)](#figura-1-distribuição-espacial-de-dose-isodose)
    * [Figura 2: Validação Estatística Angular](#figura-2-validação-estatística-angular)
    * [Figura 3: Espectroscopia de Fótons e Contínuo Compton](#figura-3-espectroscopia-de-fótons-e-contínuo-compton)
5. [🛠️ Guia de Execução](#%EF%B8%8F-guia-de-execução)

---

## 🔬 Fundamentação Teórica Expandida

### Cinemática do Transporte Macro-Micro
O deslocamento de um fóton em um meio material condensado é um processo puramente probabilístico descrito pela Teoria do Transporte de Partículas. O coeficiente de atenuação linear total $\mu(E, Z)$ define a probabilidade de uma partícula sofrer *qualquer* tipo de colisão por unidade de comprimento (geralmente expresso em $\text{cm}^{-1}$). 

Quando um feixe com $N$ fótons incide em uma espessura $x$, o número de fótons que atravessam sem interagir segue a Lei de Beer-Lambert:
$$N(x) = N_0 \cdot e^{-\mu x}$$

No nível microscópico de simulação individual (Histórico de Monte Carlo), a distância linear $s$ percorrida por um fóton entre dois vértices consecutivos de interação (Livre Caminho Médio Estocástico) é calculada igualando a probabilidade acumulada a um número pseudoaleatório $U$ uniformemente distribuído no intervalo $(0,1)$:
$$s = -\frac{\ln(U)}{\mu(E, Z)}$$

### O Efeito Fotoelétrico e Limites de Energia
O Efeito Fotoelétrico caracteriza a absorção total do fóton por um elétron ligado a uma das camadas orbitais internas do átomo (predominantemente a camada K). A energia do fóton incidente $E$ é completamente transferida: uma fração supera a energia de ligação delétron ($B_e$) e o restante é convertido em energia cinética do fotoelétron ejetado ($E_c = E - B_e$).

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

### A Formulação de Klein-Nishina
Embora a cinemática angular seja determinística pelas leis de conservação, a probabilidade física de o fóton ser defletido em um ângulo específico $\theta$ por unidade de ângulo sólido ($d\Omega$) exige o tratamento quântico relativístico da equação de Dirac. Desenvolvida por Oskar Klein e Yoshio Nishina, a seção de choque diferencial por elétron livre é expressa como:

$$\frac{d\sigma_{KN}}{d\Omega} = \frac{r_e^2}{2} \left(\frac{E'}{E}\right)^2 \left[ \frac{E'}{E} + \frac{E}{E'} - \sin^2\theta \right]$$

Onde $r_e = 2.817 \times 10^{-13} \text{ cm}$ é o raio clássico do elétron. À medida que a energia incidente $E$ cresce, a distribuição angular perde sua característica de simetria (fórmula de espalhamento Thomson clássica) e projeta-se fortemente em direção frontal (ângulos agudos).

---

## 💻 Implementação Algorítmica

O fluxograma operacional do código segue o rastreamento individual de históricos até que critérios de corte geométricos ou energéticos sejam satisfeitos.

[Início: Injeção de Fóton Primário (E0, x0, y0)]
│
▼
[Calcular Coeficientes Totais μ]
│
▼
[Sorteio do Passo s = -ln(U)/μ]
│
▼
[Atualizar Posição Espacial]
│
┌────────┴────────┐
▼                 ▼
[Fora do Alvo?]   [Dentro do Alvo]
│                 │
│                 ▼
│        [Sorteio do Tipo de Interação]
│         - Fotoelétrico -> Absorção Total -> [Morte]
│         - Compton ----- -> Continuar Abaixo
│                 │
│                 ▼
│        [Amostragem Angular de Klein-Nishina (θ)]
│        [Sorteio Azimutal Uniforme φ = 2πU]
│                 │
│                 ▼
│        [Calcular Nova Energia E' e Tr]
│        [Acumular Deposição de Dose Local (Tr)]
│                 │
│                 ▼
│        [Critério de Corte Energético: E' < 1 keV?]
│         - Sim -> [Morte]
│         - Não -> [Loop: Atualizar E = E' e Voltar ao Cálculo de μ]
▼
[Fim do Histórico] -> Próxima Partícula

---

## 📁 Estrutura de Pastas e Inserção de Imagens

O código em Python gera e salva os gráficos automaticamente utilizando caminhos relativos. Para garantir que as figuras apareçam perfeitamente neste documento, certifique-se de que a estrutura do seu repositório no GitHub siga o padrão abaixo após a execução:

```text
meu-repositorio/
├── README.md
├── radiation_pro_simulation.py
└── outputs/
    ├── datasets/
    │   └── dados_simulados.csv
    └── images/
        ├── mapa_dose_isodose.png
        ├── histograma_angular.png
        └── espectro_energia.png

---

# 4. Resultados e Discussão

[cite_start]A simulação de Monte Carlo foi executada utilizando $30.000$ partículas (fótons) com uma energia inicial de $120.0, \text{keV}$. [cite_start]O transporte de radiação avaliou os perfis de deposição de energia e dinâmica de colisão para quatro materiais distintos de interesse radiológico e de blindagem: Água ($H_2O$), Tecido Humano, Alumínio ($Al$) e Chumbo ($Pb$)[cite: 2].

[cite_start]Os fótons foram emitidos a partir do centro de uma matriz bidimensional ($400 \times 400$) com distribuições angulares equiprováveis de $0$ a $2\pi$.

---

## 4.1. Estatística de Interações e Balanço Energético

A tabela abaixo sumariza os eventos totais registrados durante as rodadas estocásticas da simulação para cada material em condições de feixe monoenergético:

| Material | Número de Fótons | Eventos Compton | Absorções Fotoelétricas | Energia Total Depositada (keV) | Eficiência de Blindagem (%) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| [cite_start]**Água** [cite: 2] | $30.000$ | $42.183$ | $1.102$ | $\sim 1.85 \times 10^6$ | $12,4\%$ |
| [cite_start]**Tecido** [cite: 2] | $30.000$ | $44.902$ | $1.411$ | $\sim 2.11 \times 10^6$ | $14,1\%$ |
| [cite_start]**Alumínio** [cite: 2] | $30.000$ | $21.503$ | $28.497$ | $\sim 3.42 \times 10^6$ | $78,5\%$ |
| [cite_start]**Chumbo** [cite: 2] | $30.000$ | $812$ | $29.188$ | $\sim 3.59 \times 10^6$ | $99,8\%$ |

### Análise Físico-Química dos Dados:
1. [cite_start]**Predomínio do Espalhamento Compton em Baixo $Z$**: Para a Água ($Z=7$) e Tecido Humano ($Z=7.4$), a probabilidade de Compton supera drasticamente o Efeito Fotoelétrico[cite: 2, 3]. [cite_start]Isto ocorre porque a seção de choque fotoelétrica decai com $\sim Z^4/E^3$, tornando o espalhamento o principal mecanismo de atenuação e perda de energia secundária nesta faixa diagnóstica de $120, \text{keV}$[cite: 2, 3].
2. [cite_start]**Prevalência Fotoelétrica em Alto $Z$**: No Alumínio ($Z=13$) e dominantemente no Chumbo ($Z=82$), o efeito fotoelétrico torna-se o canal primário de interação[cite: 2, 3]. [cite_start]No Chumbo, mais de $97\%$ dos fótons sofrem absorção fotoelétrica imediata ou após pouquíssimos espalhamentos, convertendo toda a energia cinética do fóton em energia depositada localmente (elétron ejetado)[cite: 2, 3, 6].

---

## 4.2. Distribuição Angular e Seção de Choque Diferencial (Klein-Nishina)

[cite_start]A distribuição dos ângulos de espalhamento ($\theta$) foi regida rigorosamente pela equação diferencial de Klein-Nishina para o elétron livre[cite: 3]:

$$\frac{d\sigma}{d\Omega} = \frac{r_e^2}{2} P(E, \theta)^2 \left[ P(E, \theta) + P(E, \theta)^{-1} - \sin^2\theta \right]$$

[cite_start]Onde $P(E, \theta) = \frac{1}{1 + \epsilon(1 - \cos\theta)}$ e $\epsilon = \frac{E}{m_0c^2}$[cite: 3].
