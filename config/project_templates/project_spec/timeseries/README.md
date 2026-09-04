# Timeseries Input Template

This directory is the canonical timeseries root. Copy it or generate an equivalent
layout with your own preprocessing workflow. There is no intermediate `example`
directory: each scenario is directly below this directory.

The CSV files are one-row sample templates. Replace the sample value columns and rows
with the exact zone, generator, and reserve-service names from your system. Provide
hourly, preprocessed data for every scenario and simulation year listed in
`../scenarios.csv`.

```text
{scenario}/sim_year_{n}/
├── Load/
│   ├── DAY_AHEAD_regional_Load.csv
│   └── REAL_TIME_regional_Load.csv
├── WIND/
│   ├── DAY_AHEAD_wind.csv
│   └── REAL_TIME_wind.csv
├── PV/
│   ├── DAY_AHEAD_pv.csv
│   └── REAL_TIME_pv.csv
└── Reserves/
    ├── DAY_AHEAD_regional_{product}.csv
    └── REAL_TIME_regional_{product}.csv
```

Required format:

- Header fields: `Year,Month,Day,Period` followed by value columns.
- Load value columns match `zones.csv` `load_column` values.
- Wind and PV value columns match PSY renewable component names. If a caller uses
    different names, it supplies the component-to-column mapping directly to the builder.
- Reserve value columns match PSY service names. If a caller uses different names, it
    supplies the service-to-column mapping directly to the builder.
- Values are numeric, finite, hourly, and already in the units expected by the PSY scaling-factor configuration.
- The builder wraps profiles cyclically when a forecast window crosses the end of the supplied profile.

The actual data values, number of rows, scenario directories, and reserve-product files are system-specific and must be supplied by the user.
