# Pairwise VECM results: four modes vs inventories

Logged monthly series, June 2009–December 2024 (n = 187).  
Specification: Johansen `ca.jo` with `ecdet = "trend"`, `spec = "longrun"`, `season = 12`.  
Lag K chosen by AIC over 1–12 (`ca.jo` requires K ≥ 2). VECM estimated at **r = 1** with `cajorls`.  
Significance: \* 10%, \*\* 5%, \*\*\* 1%.

R code: `estimate_pairwise_vecm.R`  
Workbook: `VECM_Results_Summary.xlsx`

---

## Specification notes

1. **Lag choice.** AIC search is 1–12, matching `load_transportation_data.R`. AIC and SBC disagree in every model (SBC prefers a much shorter VAR: 2, 2, 2, 3). AIC is used. Expanding the grid from 10 to 12 did not change K used (6, 5, 4, 5) or any α / half-life. Model 2 SBC moved from 1 to 2 because `VARselect` compares lags on a common sample of T − lag.max.
2. **Rank vs VECM.** All four VECMs are estimated at r = 1. That is the cointegrating VECM of interest, but it is only cleanly supported by **both** trace and λmax in Model 4. Models 1–2 have mixed or full-rank Johansen results; treat their long-run vectors as tentative.
3. **Half-life.** HL = ln(0.5) / ln(1 + α), reported only for significant α with −1 < α < 0. A positive α can still be error-correcting depending on the sign of β; the formula is then undefined (Model 3, inventories).
4. **Asymmetry.** YES if exactly one loading is significant at 5%, or both are significant and |α| differs by a factor of two or more.

---

## Table 1. Johansen cointegration tests (5%)

Trace 5% CV: r = 0 is 25.32; r ≤ 1 is 12.25.  
λmax 5% CV: r = 0 is 18.96; r ≤ 1 is 12.25.

| Model | Pair | K | Trace r=0 | Trace r≤1 | Trace rank | λmax r=0 | λmax r≤1 | λmax rank |
|---|---|---:|---:|---:|---|---:|---:|---|
| 1 | LTL + IC | 6 | 37.48 | 14.24 | r = 2 | 23.24 | 14.24 | r = 2 |
| 2 | Truckload + IC | 5 | 33.41 | 14.88 | r = 2 | 18.53 | 14.88 | r = 0 |
| 3 | Local + IC | 4 | 26.73 | 12.05 | r = 1 | 14.68 | 12.05 | r = 0 |
| 4 | Scheduled airfreight + IC | 5 | 32.58 | 8.19 | r = 1 | 24.40 | 8.19 | r = 1 |

**Conclusions**

- **Model 1.** Both tests reject r ≤ 1. Full rank: the two logs behave as if they are already stationary as a system. A rank-1 VECM is not supported.
- **Model 2.** Trace says full rank; λmax fails to reject r = 0 (18.53 vs 18.96). Evidence of cointegration is mixed.
- **Model 3.** Trace r = 1; λmax r = 0. Weak / mixed cointegration.
- **Model 4.** Both tests r = 1. This is the only pair with a clean cointegrating rank.

---

## Table 2. Error-correction loadings and diagnostics

| Model | α_mode | t | α_IC | t | AIC | SBC | LB p (mode) | LB p (IC) | No AC (both) |
|---|---|---:|---|---:|---:|---:|---:|---:|---|
| 1 LTL | −0.1233\*\*\* | −4.63 | −0.0089 | −0.84 | −3410.8 | −3263.7 | 0.943 | 0.998 | YES |
| 2 Truckload | −0.0461\*\*\* | −3.68 | 0.0040 | 1.19 | −3336.7 | −3202.1 | 0.365 | 0.805 | YES |
| 3 Local | −0.0089 | −0.57 | 0.0184\*\*\* | 3.26 | −3475.7 | −3353.7 | 0.003 | 0.961 | NO |
| 4 Airfreight S | −0.0499\*\*\* | −4.56 | −0.0046 | −1.35 | −3380.3 | −3245.7 | 0.642 | 0.897 | YES |

Multivariate Portmanteau (12 lags): Model 1 p = 0.64; Model 2 p = 0.47; Model 3 p = 0.017; Model 4 p = 0.077.

Inventories are **weakly exogenous** in Models 1, 2, and 4 (α_IC insignificant). Local trucking is the reverse: the **mode** is weakly exogenous and inventories adjust.

---

## Table 3. Cointegrating vectors (log mode = 1)

Relation: `log_Mode + β_IC log_IC + β_trend trend = 0`  
so `log_Mode = −β_IC log_IC − β_trend trend`. Long-run elasticity of the mode with respect to inventories is −β_IC.

| Model | β_Mode | β_IC | β_trend | LR elasticity (mode wrt IC) |
|---|---:|---:|---:|---:|
| 1 LTL | 1 | 0.128 | −0.0042 | **−0.128** |
| 2 Truckload | 1 | 1.722 | −0.0089 | **−1.722** |
| 3 Local | 1 | −1.193 | 0.0014 | **+1.193** |
| 4 Airfreight S | 1 | 1.529 | −0.0053 | **−1.529** |

**Interpretation.** In LTL, truckload, and scheduled airfreight, a long-run rise in inventories is associated with **lower** mode prices (negative elasticity). Local trucking is the opposite: inventories and local rates move together. LTL’s elasticity is small (−0.13); truckload and airfreight are above one in absolute value.

---

## Table 4. Half-lives (significant α only)

HL = ln(0.5) / ln(1 + α), in months.

| Model | Variable | α | t | Half-life |
|---|---|---|---:|---|
| 1 | log LTL | −0.1233\*\*\* | −4.63 | **5.27 months** |
| 1 | log IC | −0.0089 | −0.84 | n.a. (insignificant) |
| 2 | log truckload | −0.0461\*\*\* | −3.68 | **14.68 months** |
| 2 | log IC | 0.0040 | 1.19 | n.a. (insignificant) |
| 3 | log local | −0.0089 | −0.57 | n.a. (insignificant) |
| 3 | log IC | 0.0184\*\*\* | 3.26 | n.a. (α > 0; formula not defined) |
| 4 | log airfreight S | −0.0499\*\*\* | −4.56 | **13.55 months** |
| 4 | log IC | −0.0046 | −1.35 | n.a. (insignificant) |

LTL closes half a disequilibrium in about **five months**. Truckload and scheduled airfreight take **about 14 months**. For local trucking the significant loading is on inventories and is positive, so the AR(1) half-life formula does not apply; the implied ECT persistence 1 + α′β is about 0.97 (half-life ~22 months if that mapping is used).

---

## Table 5. Summary comparison

| Model | Trace r | λmax r | LR elast. | Faster \|α\| | Asymmetric | HL (mode) | LB clean |
|---|---|---|---:|---|---|---|---|
| 1 LTL | 2 | 2 | −0.13 | LTL rates | **YES** | 5.3 mo | YES |
| 2 Truckload | 2 | 0 | −1.72 | Truckload rates | **YES** | 14.7 mo | YES |
| 3 Local | 1 | 0 | +1.19 | Inventories | **YES** | n.a. | NO |
| 4 Airfreight S | 1 | 1 | −1.53 | Airfreight rates | **YES** | 13.5 mo | YES |

---

## Asymmetry

All four systems are **asymmetric**.

| Pattern | Models |
|---|---|
| Mode prices adjust; inventories weakly exogenous | 1 LTL, 2 truckload, 4 scheduled airfreight |
| Inventories adjust; mode weakly exogenous | 3 local trucking |

That is the economically useful split. Line-haul and airfreight rates do the error-correction; local trucking does not.

---

## Model-by-model findings

**Model 1 — LTL.** Fastest estimated adjustment (α = −0.12\*\*\*, 5.3 months) and clean residuals, but Johansen says **r = 2**. Do not treat this as a standard cointegrated VECM. The small long-run elasticity (−0.13) is also a warning that the vector is weakly identified.

**Model 2 — Truckload.** Mode-driven adjustment (14.7 months), inventories exogenous, residuals clean. Rank tests disagree (trace r = 2, λmax r = 0). Long-run elasticity −1.72: a 1% inventory increase is associated with a 1.7% decline in truckload rates in the long run, if the rank-1 relation is maintained.

**Model 3 — Local trucking.** Only mixed evidence of cointegration, and the local-rate equation does not error-correct. Inventories adjust (α_IC = +0.018\*\*\*). Ljung-Box rejects white residuals in the local equation. Weakest specification of the four.

**Model 4 — Scheduled airfreight.** Best-behaved pair: both rank tests give r = 1, airfreight α = −0.050\*\*\* (13.5 months), inventories weakly exogenous, Ljung-Box clean. Long-run elasticity −1.53. Preferred VECM for further work (IRFs, VECM forecasting).

---

## Short-run lags (1–3), selected significant terms

- **LTL:** own lag 1 ΔLTL = −0.250\*\*\*; ΔIC lag 1 in the LTL equation = 0.623\*\*\*. IC equation: own lags 1 and 3 significant.
- **Truckload:** ΔIC lag 1 in the truckload equation = 0.822\*\*\*. IC equation: truckload lags 1 and 3 significant.
- **Local:** Δlocal lag 3 in the local equation = 0.183\*\*. IC equation: own lags 1 and 3, plus local lag 2.
- **Airfreight S:** own lags 1 (−0.136\*) and 3 (−0.173\*\*). IC equation: airfreight lag 2 = 0.088\*\*\* and own lags 1 and 3.

Full coefficient lists are in `VECM_Results_Summary.xlsx` (sheet `6_VECM_all_coefs` and `lags_1to3`) and `vecm_output/`.

---

## Re-run after lag.max = 12 and table-join fix

`estimate_pairwise_vecm.R` now searches AIC over 1–12 and joins diagnostics on `model`/`pair` before printing Table 2 (so AIC/SBC are not recycled).

Primary K, α, and half-lives are unchanged:

| Model | K | α_mode | α_IC | HL (mode) | AIC | SBC |
|---|---:|---|---|---|---:|---:|
| 1 LTL | 6 | −0.1233*** | −0.0089 | 5.27 | −3410.79 | −3263.66 |
| 2 Truckload | 5 | −0.0461*** | 0.0040 | 14.68 | −3336.65 | −3202.08 |
| 3 Local | 4 | −0.0089 | 0.0184*** | n.a. | −3475.65 | −3353.69 |
| 4 Airfreight S | 5 | −0.0499*** | −0.0046 | 13.55 | −3380.29 | −3245.72 |

AIC/SBC differ by row (join fix confirmed).

Robustness (`estimate_robustness_vecm.R`, same 1–12 grid): some warehouse pairs were previously capped at K = 10. New AIC K is 11 for local vs both proxies and 12 for truckload vs storage. Truckload vs storage mode α falls from −0.076*** to −0.038* and the Johansen rank for that pair drops to r = 0. Local IC-side adjustment still holds (α_IC-proxy = 0.018***/0.072***/0.048**). Verdict unchanged: LTL, truckload, and airfreight are not robust; local is partially robust.

Details: `outputs/robustness/ROBUSTNESS_SUMMARY.txt`.
