# Transportation VECM

Monthly U.S. transportation cost indexes and inventory levels (2009–2024) for a vector error-correction model (VECM) of transportation costs and inventory.

## Data

`Transportation_Inventory_Complete_with_WPU30.xlsx` (Jun 2009–Dec 2024, 187 months, no missing values)

- **IC** = `Total_Inventories` (Census MTIS, $ millions)
- **TC** = `TC_Aggregate_WPU30` (BLS WPU30 transportation services aggregate)
- **Modes** = seven BLS transportation service cost indexes, including `Warehousing_Storage`

## Analysis

`load_transportation_data.R` uses natural logs (not raw levels):

1. ADF tests (`urca::ur.df`) on log levels and log differences
2. Johansen cointegration (`urca::ca.jo`) of log IC vs log TC and vs each mode

```r
Rscript load_transportation_data.R
```

Plots and Johansen tables are written to `plots/`.
