# Transportation VECM

Monthly U.S. transportation cost indexes and inventory levels (2003–2024) for a vector error-correction model (VECM) of transportation costs and inventory.

## Data

`Final_Transportation_Inventory_Dataset .xlsx`

- **IC** = `Total_Inventories` (Census MTIS, $ millions)
- **TC modes** = seven BLS transportation service cost indexes, including `Warehousing_Storage`

## Analysis

`load_transportation_data.R` loads the Excel file, plots levels, and runs Augmented Dickey-Fuller tests on levels and first differences.

```r
Rscript load_transportation_data.R
```

Plots are written to `plots/`.
