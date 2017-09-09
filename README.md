# Robust Monetary Policy under Time-Varying Model Uncertainty

This repo contains the code used for the thesis <a href="http://nbviewer.jupyter.org/github/christophermosch/robustness/blob/master/Robust%20Monetary%20Policy%20under%20Time-Varying%20Model%20Uncertainty.pdf" target="_blank">Robust Monetary Policy under Time-Varying Model Uncertainty</a>. Below is a brief description for replicating the results.

## Simulations in Section 5

`Sim_.m` simulates the evolution of the New Keynesian model used in Section 5. It allows for various settings. First, one can choose between constant and time-varying uncertainty on the side of the central bank (theta\_const = 0 or 1). Moreover, one can either explore the long-run behavior of the model (IRF = 0) or the impulse responses (IRF = 1 or 2). In the former case, a normally i.i.d. shock occurs every period. When impulse responses are selected, a shock of one standard deviation occurs in the first period. Finally, the code includes the option to compute the detection error probabilities for the current approximating model, loss function, discount factor, and sample size, i.e. for the values chosen for α, β, γ, ρ<sub>1</sub>, ρ<sub>2</sub>, σ<sub>1</sub>, σ<sub>2</sub>, λ<sub>y</sub>, λ<sub>i</sub>, and t. Figures 3 and A.1 were created with `Sim_.m`.

`Sim_IRF.m` repeatedly simulates the impulse responses under time-varying uncertainty in order to obtain their distribution. It was used to create Figures 4 and 5.

`Breakdown.m` calculates the breakdown point given the approximating model, the planner’s loss function, and his discount factor. It does so by minimizing θ under the condition that the matrix V<sub>t</sub> in the value function always has to be positive-definite while iterating backwards to convergence.

`BreakLoss.m` is the function that is minimized by fmincon in `Breakdown.m`.

`DetErrProb.m` is called by `Sim_m` to calculate the detection error probabilities. The function is as in Giordani and Söderlind (2004) , except that it uses θ instead of -1/θ.

`Var1SimPs.m` is called by `DetErrProb.m` to model the evolution of the state variables. The function is as in Giordani and Söderlind (2004). 

`DiscAlgR.m`. The function transforms the decision problem such that the rational expectations techniques of Backus and Driffil (1986) can be applied to iterate backwards to convergence. It calls on `DiscAlg.m`. The function is as in Giordani and Söderlind (2004), except that it additionally checks whether the matrix V<sub>t</sub> is positive-definite during each iteration.

`DiscAlg.m`. Iterates backwards until convergence. To do so, it calls `DiscAlg2.m` during each iteration. The function is as in Giordani and Söderlind (2004), except that it additionally checks whether the matrix V<sub>t</sub> is positive-definite during each iteration. 

`DiscAlg2.m`. Calculates the left-hand sides of equations (3.9) – (3.12). The function is as in Giordani & Söderlind (2004).

## Inference in Section 6

`Emp_IRF_hist.m` Estimates the impulse responses for the rolling samples created from the input data. To this end, it calls on impulseresponse.m. Figure 6 was created with `Emp_IRF_hist.m`.

`Comparison_IRF_hist.m` combines `Sim_IRF.m` and `Emp_IRF_hist.m` to directly compare the predicted distribution of the impulse responses with the empirical distribution. The code was used to create Figures A.1 and A.2.

`Emp_Var_order_selection.m` is used to select the number of lags in the VAR from which the impulse responses are computed. It was created by Kevin Sheppard and is available [here](https://www.kevinsheppard.com/MFE\_MATLAB#Vector\_Autoregressions).

`impulseresponse.m`, `vectorar.m`, `newlagmatrix.m`, `vectorarvcv.m` are called to estimate a VAR and the impulse responses. They were created by Kevin Sheppard and are part of his [toolbox](https://www.kevinsheppard.com/MFE\_Toolbox).

## References
- Backus, D., & Driffill, J. (1986). The consistency of optimal policy in stochastic rational expectations models (CEPR Discussion Paper No. 124). C.E.P.R. Discussion Papers.
- Giordani, P., & Söderlind, P. (2004). Solution of macromodels with Hansen–Sargent robust policies: some extensions. Journal of Economic Dynamics and Control, 28(12), 2367–2397.