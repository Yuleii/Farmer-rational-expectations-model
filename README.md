
This repository is aimed to solve the multivariate rational expectations model correspond to pages 41 - 51 of the book
by [Roger E. A. Farmer, "Macroeconomics of Self-fulfilling Prophecies", 2nd edition, MIT Press, 1999](https://mitpress.mit.edu/books/macroeconomics-self-fulfilling-prophecies-second-edition).

### Equilibrium conditions of the stochastic Ramsey economy

<img src="https://render.githubusercontent.com/render/math?math=\frac{1}{c_{t}}=\beta E_{t}\left[\frac{1}{c_{t%2B1}}\left\{(1-\delta)%2BF_{k}^{t%2B1}\right\}\right]">

<img src="https://render.githubusercontent.com/render/math?math=k_{t%2B1}=(1-\delta) k_{t}%2Bs_{t} k_{t}^{\alpha}-c_{t}">

<img src="https://render.githubusercontent.com/render/math?math=s_{t%2B1}=s_t^{\rho}v_{t%2B1}">

where: <img src="https://render.githubusercontent.com/render/math?math=F_k^{t%2B1}=\alpha s_{t%2B1} k_{t%2B1}^{\alpha-1}">


### Disturbances

log of disturbances to productivity is normally distributed <img src="https://render.githubusercontent.com/render/math?math=\log v_{t%2B1} \sim \mathcal{N}(0,\sigma_v^2)">, with <img src="https://render.githubusercontent.com/render/math?math=\sigma_v^2=0.007^2">

### Parameters
<img src="https://render.githubusercontent.com/render/math?math=\alpha = 0.3">

<img src="https://render.githubusercontent.com/render/math?math=\beta = 0.95">

<img src="https://render.githubusercontent.com/render/math?math=\delta = 0.1">

<img src="https://render.githubusercontent.com/render/math?math=\rho = 0.9">


Given inital 10 periods:

<img src="https://render.githubusercontent.com/render/math?math=\{ \log s_t\}_{t=0}^{9}=\{0, -0.005, -0.009, -0.013, -0.022, -0.021, -0.019, -0.011, -0.012, -0.003\}">




