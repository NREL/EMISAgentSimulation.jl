# Project Initialization Input Set

This directory is the canonical project-init input contract. Copy it, edit the CSV
values and paths, and pass the copied directory to the project initializer. The
initializer expects these filenames exactly; users may also produce the same files from
an independent preprocessing workflow.

The filenames in this directory are the initializer's input filenames. Do not rename them.

Required tables:

- `project_spec.csv`
- `zones.csv`
- `technologies.csv`
- `scenarios.csv`
- `projectoptions.csv`
- `system_config.csv`

`projectexisting.csv`, `devices_to_remove.csv`, and `outages.csv` are optional. The PSY system JSON and actual profile CSVs are external files referenced by these inputs and are not represented by placeholder data here. `reference_case_dir` is optional and may remain blank for a first-time project.

The `project_spec.csv` outage path may point to `outages.csv` or to a supported long-form or wide outage dataset. `start_year` and `pcm_scenario` set the defaults rendered into each case template; they may be changed for an individual case with `--override`. Timeseries are discovered from the canonical folder and filename layout below; users are responsible for placing their preprocessed files there.

## External input layout

The paths in `project_spec.csv` may be absolute or relative to the spec directory.
The referenced system and profile data should have this structure:

```text
external_inputs/
├── system.json
├── outages.csv                         # or a supported long-form/wide outage file
└── timeseries/
	└── {scenario}/sim_year_{n}/
		├── Load/
		│   ├── DAY_AHEAD_regional_Load.csv
		│   └── REAL_TIME_regional_Load.csv
		├── WIND/
		│   ├── DAY_AHEAD_wind.csv
		│   ├── REAL_TIME_wind.csv
		├── PV/
		│   ├── DAY_AHEAD_pv.csv
		│   └── REAL_TIME_pv.csv
		└── Reserves/
			├── DAY_AHEAD_regional_{product}.csv
			└── REAL_TIME_regional_{product}.csv
```

Each profile CSV must contain a header row and numeric component or service columns.
The canonical profile contract uses hourly, preprocessed data. Load columns correspond
to `zones.csv` `load_column` values; wind and PV columns correspond to PSY component
names; reserve columns correspond to PSY service names. The first four timestamp fields
may be retained as `Year,Month,Day,Period`; the builder uses the remaining value columns.
If source column names differ from PSY names, provide the corresponding column mapping
when calling the builder. Profile preprocessing, unit conversion, and weather-year
creation remain the user's responsibility.

The minimal outage template uses `GEN_UID,FOR,MTTR Hr`. Existing projects may instead
reference the supported outage formats used by the extractor, provided their path is
set in `project_spec.csv`.

The package-level path defaults are maintained in `config/timeseries_defaults.csv` with
the columns `kind,market_stage,directory,filename`. The builder uses these defaults
when a profile path is not explicitly supplied. A project or caller may provide an
alternate copy through `timeseries_defaults_file`; explicit profile paths still take
priority. This file is configuration, not a second required user input.

## Reference case

`reference_case_dir` is a convenience source for compatible market, investor, finance,
and other downstream files. It is not required. When it is blank, initialization must
use the repository templates, `config/project_defaults.csv`, and deterministic defaults
where applicable. A new project still must provide its own PSY system, zones,
technologies, scenarios, project options, and timeseries inputs.
