# HDF5 Serialization Plan

**Goal:** Replace all `FileIO.save` / `FileIO.load` / `JLD2.jldsave` calls in EMIS with
field-by-field HDF5 (`.h5`) saves. `PSY.System` and PRAS systems continue to use
`PSY.to_json` (unchanged).

**Implementation file:** `src/utils/save_load_utils.jl` (created; not yet wired into call sites)

---

## 1. Affected call sites

| File | Lines | What is saved | New filename |
|------|-------|---------------|--------------|
| `investment_simulation.jl` | 142–149, 471–478 | `derating_factors` per scenario/year | `derating_data_year_N.h5` |
| `investment_simulation.jl` | 375–382 | `shortfall_data` per scenario/year | `shortfall_data_<scenario>_year_N.h5` |
| `investment_simulation.jl` | 500–506 | `simulation_data` per year | `simulation_data_year_N.h5` |
| `investment_simulation.jl` | 507–513 | `clean_energy_percentage` per year | `clean_energy_percentage_year_N.h5` |
| `investment_simulation.jl` | 515–521 | Duplicate JLD2 saves of the above | **Delete — redundant** |
| `investment_simulation.jl` | 556–560 | Final `clean_energy_percentage` | `clean_energy_percentage.h5` |
| `investment_simulation.jl` | 561–565 | Final `simulation_data` | `simulation_data.h5` |
| `ra_utils.jl` | 289–295 | 7 debug scalars (temp data, marked TODO) | Consolidate → `temp_debug_data.h5` |
| `actual_market_simulation.jl` | 313–344 | Realized market data (~20 fields) | `year_N.h5` |
| `expected_market_simulation.jl` | 113–166 | Expected market data (~30 fields), investor dir | `<scenario>_year_N.h5` |
| `expected_market_simulation.jl` | 168–182 | Same data, sys_results dir | `<scenario>_year_N.h5` |
| `investor_iteration.jl` | 40–41 | **LOAD** of expected market data | Read from `<scenario>_year_N.h5` |

PSY systems are already saved correctly with `PSY.to_json` at `investment_simulation.jl:523–533`
and `ra_utils.jl:296` — these are **not changed**.

**Dependency changes (Project.toml):**
- Add: `HDF5`
- Remove: `FileIO`, `JLD2`

---

## 2. Global encoding rules

Every Julia type maps to HDF5 as follows. These rules are applied uniformly throughout
`save_load_utils.jl`.

| Julia type | HDF5 encoding |
|-----------|---------------|
| `String` | Scalar string dataset |
| `Int64`, `Float64`, `Bool` | Scalar numeric dataset |
| `Vector{T}` | 1-D dataset |
| `Matrix{T}` | 2-D dataset |
| `Symbol` | String dataset — `string(sym)` on save, `Symbol(s)` on load |
| `NamedTuple{(:a,:b),...}` | Group with datasets named `"a"`, `"b"` |
| `Dict{String, V}` | Group; one child group/dataset per key |
| `Dict{Symbol, V}` | Group with parallel arrays `keys::String[]` + `values::V[]` |
| `Dict{Tuple{String,String}, Vector{Float64}}` | Group with `key_a::String[]`, `key_b::String[]`, `values/<i>::Float64[]` |
| `AxisArray{T,N}` | Group: `data` (N-D array) + `axis<i>_name`, `axis<i>_values`, `axis<i>_type` per axis |
| `Union{Nothing, T}` | HDF5 attribute `is_nothing = true` on group if `Nothing`; otherwise write value normally |
| Abstract type (BuildPhase, Project, Forecast, RiskPreference, Product) | HDF5 attribute `"type" = "ConcreteTypeName"` on group — used to dispatch on load |
| `DataFrame` | Group: `col_names::String[]` + one dataset per column |
| `PSY.System` | **Excluded from .h5** — use `PSY.to_json` companion file |
| `PSY.*GenerationCost` | Group with extracted scalar fields + cost-curve breakpoint arrays |

A root attribute `schema_version = 1` is written to every file for future migration support.

---

## 3. `simulation_data_year_N.h5` — full layout

This is the main checkpoint file written once per investor iteration year. It stores the
entire `AgentSimulation` struct except PSY systems.

```
simulation_data_year_N.h5
│
├── [attr] schema_version = 1
│
├── results_dir              String
├── iteration_year           Int64
├── rep_period_interval      Int64
├── carbon_tax               Float64[N_years]
├── rec_requirement          Float64[N_years]
├── zones                    String[N_zones]
│
├── markets/                 Dict{Symbol, Bool}
│   ├── keys                 String[N]
│   └── values               Bool[N]
│
├── lines/                   Vector{ZonalLine} — stored as parallel arrays
│   ├── names                String[N]
│   ├── from_zone            String[N]
│   ├── to_zone              String[N]
│   └── active_power_limit   Float64[N]
│
├── rep_periods/             Dict{String, Dict{Int64, Dict{Int64,Int64}}}
│   └── <scenario>/
│       └── <year>/
│           ├── keys         Int64[N]
│           └── values       Int64[N]
│
├── hour_weight/             Dict{String, Dict{Int64, Vector{Float64}}}
│   └── <scenario>/
│       └── <year>           Float64[N_hours]
│
├── peak_load/               Dict{String, Dict{Int64, Float64}}
│   └── <scenario>/
│       ├── keys             Int64[N_years]
│       └── values           Float64[N_years]
│
├── derating_data/           Dict{String, DataFrame}
│   └── <tech_type>/
│       ├── col_names        String[N_cols]
│       └── <col_name>       typed array per column
│
├── resource_adequacy/       Dict{String, ResourceAdequacy}
│   └── <scenario>/
│       ├── target_keys      String[N]
│       ├── target_values    Float64[N]
│       ├── delta_irm        Float64[N_years]
│       └── metrics/
│           └── <year>/
│               ├── keys     String[N]
│               └── values   Float64[N]
│
└── investors/               Vector{Investor}
    └── <investor_name>/
        ├── name                         String
        ├── data_dir                     String
        ├── avg_block_size               Int64
        ├── fixed_block_size             Bool
        ├── rep_period_interval          Int64
        ├── cap_cost_multiplier          Float64
        ├── max_annual_projects          Int64
        ├── retirement_lookback          Int64
        ├── carbon_tax                   Float64[N_years]
        ├── markets                      String[N]   (Symbol → String)
        │
        ├── preference_multiplier_range/
        │   ├── min                      Float64
        │   └── max                      Float64
        │
        ├── risk_preference/
        │   ├── [attr] type = "RiskNeutral" | "RiskAverse"
        │   └── (RiskAverse only:)
        │       ├── constant             Float64
        │       ├── multiplier           Float64
        │       └── risk_coefficient     Float64
        │
        ├── forecast/
        │   ├── [attr] type = "Perfect" | "Imperfect"
        │   ├── [attr] kalman_filter_is_nothing = true|false  (Imperfect only)
        │   ├── kalman_filter/           (if Imperfect and not nothing)
        │   │   ├── investor_belief/
        │   │   │   ├── process_covariance/      AxisArray{Float64,1}
        │   │   │   │   ├── data                 Float64[N]
        │   │   │   │   ├── axis1_name           String
        │   │   │   │   ├── axis1_values         values array
        │   │   │   │   └── axis1_type           String
        │   │   │   ├── measurement_covariance/  same layout
        │   │   │   ├── state_transition_matrix  Float64[N,N]
        │   │   │   └── state_measurement_matrix Float64[M,N]
        │   │   ├── state_estimate/              AxisArray{Float64,1}
        │   │   └── error_covariance_estimate/   AxisArray{Float64,2}
        │   └── scenarios/
        │       └── <scenario_name>/
        │           ├── name                     String
        │           ├── probability              Float64
        │           ├── parameter_multipliers/
        │           │   ├── [attr] is_nothing    Bool
        │           │   ├── keys                 String[N]
        │           │   └── values               Float64[N]
        │           └── parameter_values/
        │               └── <year>/              AxisArray{Float64,2}
        │                   ├── data             Float64[N_products, N_years]
        │                   ├── axis1_name/values/type
        │                   └── axis2_name/values/type
        │
        ├── rep_hour_weight/             Dict{String, Dict{Int64, Vector{Float64}}}
        │   └── <scenario>/<year>        Float64[N_hours]
        │
        ├── chron_weights/               Dict{String, Dict{Int64, Matrix{Int64}}}
        │   └── <scenario>/<year>        Int64[N_rows, N_cols]
        │
        ├── portfolio_preference_multipliers/
        │   ├── key_a                    String[N_entries]
        │   ├── key_b                    String[N_entries]
        │   └── values/
        │       └── <flat_index>         Float64[N_years]
        │
        ├── market_prices/
        │   ├── energy_price/            Dict{String, AxisArray{Float64,3}}
        │   │   ├── [attr] is_nothing    Bool
        │   │   └── <scenario>/          AxisArray group (data + 3 axes)
        │   ├── reserve_price/           Dict{String, Dict{String, Array{Float64,2}}}
        │   │   ├── [attr] is_nothing    Bool
        │   │   └── <product>/<scenario> Float64[N_hours, N_years]
        │   ├── capacity_price/          Dict{String, AxisArray{Float64,1}}
        │   │   ├── [attr] is_nothing    Bool
        │   │   └── <scenario>/          AxisArray group (data + 1 axis)
        │   ├── rec_price/               same layout as capacity_price
        │   └── inertia_price/           Dict{String, AxisArray{Float64,2}}
        │       ├── [attr] is_nothing    Bool
        │       └── <scenario>/          AxisArray group (data + 2 axes)
        │
        └── projects/
            └── <project_name>/
                ├── [attr] project_type  "ThermalGenEMIS"|"BatteryEMIS"|"RenewableGenEMIS"|"HydroGenEMIS"
                ├── [attr] build_phase   "Existing"|"Planned"|"Queue"|"Option"|"Retired"
                ├── name                 String
                ├── decision_year        Int64
                ├── construction_year    Int64
                ├── retirement_year      Int64
                ├── end_life_year        Int64
                │
                ├── tech/
                │   ├── [attr] tech_type  "ThermalTech"|"BatteryTech"|"RenewableTech"|"HydroTech"
                │   ├── type              String
                │   ├── bus               String
                │   ├── zone              String
                │   ├── FOR               Float64
                │   ├── MTTR              Int64
                │   ├── active_power_limits/    min, max  Float64
                │   ├── ramp_limits/            [attr]is_nothing, up, down  Float64
                │   │
                │   ├── (ThermalTech / HydroTech only:)
                │   │   └── time_limits/        [attr]is_nothing, up, down  Float64
                │   │
                │   ├── (ThermalTech only:)
                │   │   ├── fuel                String
                │   │   ├── fuel_cost           Float64
                │   │   └── heat_rate_curve/    x::Float64[N], y::Float64[N]
                │   │
                │   ├── (BatteryTech only:)
                │   │   ├── input_active_power_limits/   min, max
                │   │   ├── output_active_power_limits/  min, max
                │   │   ├── storage_capacity/            min, max
                │   │   ├── storage_level_limits/        min, max
                │   │   ├── initial_storage_capacity_level  Float64
                │   │   ├── rating               Float64
                │   │   ├── soc                  Float64
                │   │   ├── base_power           Float64
                │   │   └── efficiency/          in, out  Float64
                │   │
                │   └── operation_cost/
                │       ├── [attr] is_nothing    Bool
                │       ├── [attr] cost_type     "ThermalGenerationCost"|"RenewableGenerationCost"|"HydroGenerationCost"
                │       ├── fixed                Float64
                │       ├── shut_down            Float64    (Thermal only)
                │       ├── variable/            [attr]is_nothing, x::Float64[N], y::Float64[N]
                │       └── start_up/            hot, warm, cold  Float64  (Thermal only)
                │
                ├── finance_data/
                │   ├── investment_cost           Float64[N_years]
                │   ├── effective_investment_cost Float64
                │   ├── preference_multiplier     Float64[N]
                │   ├── lag_time                  Int64
                │   ├── life_time                 Int64
                │   ├── capex_years               Int64
                │   ├── fixed_OM_cost             Float64
                │   ├── queue_cost                Float64[N_years]
                │   ├── discount_rate             Float64
                │   ├── expected_npv              Float64[N_years]
                │   ├── expected_utility          Float64[N_years]
                │   ├── annual_cashflow           Float64[N_years]
                │   ├── ownedby                   String
                │   ├── scenario_total_utilization/<scenario>   Float64[N_prod, N_years]
                │   ├── scenario_npv/<scenario>                 Float64[N_years]
                │   ├── scenario_utility/<scenario>             Float64[N_years]
                │   ├── scenario_profit/<scenario>/<year>/      AxisArray{Float64,2}
                │   └── realized_profit/                        AxisArray{Float64,2}
                │
                └── products/
                    └── <product_name>/           (Symbol → String as group name)
                        ├── [attr] product_type   "Energy"|"Capacity"|"OperatingReserve{...}"|"REC"|"Inertia"
                        ├── name                  String  (Symbol → String)
                        │
                        ├── (Energy:)       marginal_cost, expected_production, capacity_factors/<scen>
                        ├── (Capacity:)     capacity_bid, derating/keys+values, accepted_perc/<scen>
                        ├── (OperatingReserve:) max_limit, marginal_cost
                        ├── (REC:)          expected_certificates, correction_factor, rec_bid
                        └── (Inertia:)      synchronous, h_constant, marginal_cost
```

---

## 4. Companion JSON files (PSY systems — unchanged)

These sit alongside the `.h5` file and are already implemented:

```
simulation_data_year_N.h5          ← everything above
sys_MD_year_N.json                 ← PSY.to_json(sys_MDs[N], ...)
sys_UC_year_N.json                 ← PSY.to_json(sys_UCs[N], ...)
sys_ED_year_N.json                 ← PSY.to_json(sys_EDs[N], ...)
sys_PRAS_<scenario>_year_N.json   ← PSY.to_json(sys_PRAS[scenario], ...)
```

---

## 5. Implementation — function pairs in `save_load_utils.jl`

| # | Save function | Load function | Handles |
|---|--------------|--------------|---------|
| 1 | `save_simulation` | `load_simulation` | `AgentSimulation` — top-level entry point |
| 2 | `save_investor!` | `load_investor` | `Investor` |
| 3 | `save_project!` | `load_project` | `ThermalGenEMIS`, `BatteryEMIS`, `RenewableGenEMIS`, `HydroGenEMIS` |
| 4 | `save_tech!` | `load_tech` | `ThermalTech`, `BatteryTech`, `RenewableTech`, `HydroTech` |
| 5 | `save_finance!` | `load_finance` | `Finance` |
| 6 | `save_product!` | `load_product` | `Energy`, `Capacity`, `OperatingReserve{T}`, `REC`, `Inertia` |
| 7 | `save_market_prices!` | `load_market_prices` | `MarketPrices` |
| 8 | `save_forecast!` | `load_forecast` | `Perfect`, `Imperfect` |
| 9 | `save_kalman_filter!` | `load_kalman_filter` | `KalmanFilter`, `InvestorBelief` |
| 10 | `save_scenario!` | `load_scenario` | `Scenario` |
| 11 | `save_resource_adequacy!` / `save_resource_adequacy_dict!` | `load_resource_adequacy` / `load_resource_adequacy_dict` | `ResourceAdequacy` |
| 12 | `save_operation_cost!` | `load_operation_cost` | PSY `*GenerationCost` types |
| 13 | `save_axisarray!` | `load_axisarray` | `AxisArray{T,N}` — generic helper used by many of the above |
| 14 | `save_dataframe!` | `load_dataframe` | `DataFrames.DataFrame` — generic helper |

Helper functions used internally:

- `save_lines!` / `load_lines` — `Vector{ZonalLine}` as parallel arrays
- `save_rep_periods!` / `load_rep_periods` — nested int dict
- `save_nested_dict_vf!` / `load_nested_dict_vf` — `Dict{String, Dict{Int64, Vector{Float64}}}`
- `save_nested_dict_f!` / `load_nested_dict_f` — `Dict{String, Dict{Int64, Float64}}`
- `save_dict_symbol_bool!` / `load_dict_symbol_bool` — `Dict{Symbol, Bool}`
- `save_chron_weights!` / `load_chron_weights` — `Dict{String, Dict{Int64, Matrix{Int64}}}`
- `save_portfolio_preference_multipliers!` / `load_portfolio_preference_multipliers`
- `save_risk_preference!` / `load_risk_preference`
- `save_derating_data!` / `load_derating_data`

---

## 6. Next steps

1. **Wire call sites** — replace `FileIO.save` / `FileIO.load` / `JLD2.jldsave` in:
   - `src/investment_simulation.jl`
   - `src/resource_adequacy/ra_utils.jl`
   - `src/markets_simulation/actual_market_simulation.jl`
   - `src/markets_simulation/expected_market_simulation.jl`
   - `src/investor_functions/investor_iteration.jl`

2. **Update `Project.toml`** — add `HDF5`, remove `FileIO` and `JLD2`.

3. **Test round-trip** — save then load a real `AgentSimulation` from a completed year and
   verify field equality.

4. **Handle `CaseDefinition`** — `load_simulation` currently returns a placeholder
   `CaseDefinition()`; the caller must repopulate it. Either extend `save_load_utils.jl`
   to serialize `CaseDefinition`, or reconstruct it from the run script inputs.
