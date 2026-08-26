# Transportation VECM

Monthly U.S. transportation cost indexes and inventory levels (2009–2024) for a vector error-correction model (VECM) of transportation costs and inventory.

## Data

`Transportation_Inventory_Complete_10_Modes_WPU30.xlsx` (Jun 2009–Dec 2024, 187 months, no missing values)

- **IC** = `Total_Inventories` (Census MTIS, $ millions)
- **TC** = `TC_Aggregate_WPU30` (BLS WPU30 transportation services aggregate)
- **Modes** = ten BLS transportation service cost indexes, including `Rail_Transportation`, `Inland_Water_Freight`, `Deep_Sea_Freight`, and `Warehousing_Storage`

## Analysis

`load_transportation_data.R` uses natural logs (not raw levels):

1. ADF tests (`urca::ur.df`) on log levels and log differences
2. Johansen cointegration (`urca::ca.jo`) of log IC vs log TC and vs each mode

```r
Rscript load_transportation_data.R
```

Each run clears `plots/` and writes only the current figures and Johansen tables.

## Pairwise VECMs

`estimate_pairwise_vecm.R` estimates four rank-1 VECMs (LTL, truckload, local trucking, scheduled airfreight vs inventories). Lag K is AIC over 1–12 (`ca.jo` requires K ≥ 2).

```r
Rscript estimate_pairwise_vecm.R
```

Outputs: `VECM_Results_Summary.xlsx`, `VECM_Results_Report.md`, and `vecm_output/`.
