# SablefishMSE - Management strategy evaluation for Alaska Sablefish

**Joshua A. Zahner (_UAF_)**, Ben Williams (_NOAA_), Curry Cunningham (_UAF_), Dan Goethel (_NOAA_), Matt Cheng (_UAF_), Pete Hulson (_NOAA_), Chris Lunsford (_NOAA_) 


A management strategy evaluation simulation framework for assessing alternative management options for Alaska sablefish in the North Pacific ocean, under the jurisidiction of the U.S. North Pacific Fisheries Management Council. 

The MSE operating model (OM) is built using the [`afscOM` R package](https://github.com/BenWilliams-NOAA/afscOM), a generalized fisheries operating model implementation. The OM is an age-structured, multi-sex, single region model, with two active fisheries and two scientific surveys, built with the same demographic parameters as are used/estimated by the [2023 Alaska sablefish stock assessment](https://github.com/dgoethel/2023-Sablefish-SAFE/). The MSE estimation model (EM) is a modified version of the [`SpatialSablefishAssessment` TMB model](https://github.com/Craig44/SpatialSablefishAssessment) built in 2022, that was updated to include a recruitment bias ramp and to fit to sex-disaggregated age composition data.

A range of recruitment models and harvest control rules (HCRs) are implemented to allow for testing the efficacy of many different management strategies across a range of reasonable future recruitment scenarios.

Examples of how to run the full MSE simulation loop is available at `dev/sablefish_mse_example.r`

### Project background and objectives
---
Alaskan sablefish (_Anoplopoma fimbria_) are currently managed using the North Pacific Fishery Management Council’s (NPFMC) $F_{40}$ harvest control rule (HCR). However, sablefish are a long-lived, relatively slow growing species and generic HCRs aimed at maximizing yearly harvest (e.g., spawner-per-recruit, SPR, based maximum sustainable yield proxies) may not perform adequately for achieving key conservation and fishery performance metrics (e.g., maintaining a robust age structure and maximizing long-term fishery yield). To address scientific and stakeholder concerns regarding the robustness of the current HCR for sablefish, a closed loop simulation tool will be developed and implemented to test the efficacy and robustness of current and alternate HCRs as well as spawning metrics through management strategy evaluation (MSE; Punt et al. 2016). The aim of the study will be to identify HCRs that can achieve both conservation and economic priorities, while also exploring how assumptions regarding calculation of spawning potential impact HCR robustness.





