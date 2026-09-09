# Project Initialization Guide

Project initialization turns a user-provided system input bundle into the directory
tree EMIS needs to run cases. It is meant for a user who has a PSY/Sienna system,
preprocessed time series, outage data, and a small set of CSVs describing zones,
technologies, scenarios, investors, and build options.

The starter input bundle lives at `config/project_templates/project_spec/`. Copy that
folder, edit the CSVs, add your real system data, then run `initialize_emis_project`.

## Basic Workflow

1. Copy `config/project_templates/project_spec/` to a user input folder.
2. Edit `project_spec.csv` with your system name, paths, investors, and generated
   directory names.
3. Edit `zones.csv`, `technologies.csv`, `scenarios.csv`, and `projectoptions.csv`.
4. Add your PSY system JSON, its HDF5 time-series sidecar, outage file, and real
   time series tree.
5. Add optional overrides only when your study needs custom market, finance, or
   case settings.
6. Run `initialize_emis_project(input_dir; output_dir=project_root)`.
7. Run `new_emis_case(project_root, case_name; overrides=...)` for each scenario run.

Keep the user input folder separate from generated output. The input folder should
contain only files the user provides or edits.

## Required User Inputs

These files must be in the copied input folder:

- `project_spec.csv`: paths, investors, and generated directory names
- `zones.csv`: EMIS zones and load-column names
- `technologies.csv`: unit-type classifications and technology defaults
- `scenarios.csv`: scenario names, PCM labels, and weather years
- `projectoptions.csv`: minimal new-build candidates
- `system_config.csv`: runtime reader settings

These data files must exist either in the input folder or at paths referenced from
`project_spec.csv`:

- PSY/Sienna system JSON from `system_filepath`
- matching `<system_stem>_time_series_storage.h5` sidecar
- preprocessed time series tree from `time_series_data_dir`
- outage data from `outage_filepath`, unless the study intentionally uses defaults

Use `config/project_templates/project_spec/timeseries/README.md` for the required
time series folder layout.

## Optional User Inputs

Add these only when needed:

- `projectexisting.csv`: existing-project ownership. Omit it for all-new-entrant
  cases or for investors with no existing projects.
- `devices_to_remove.csv`: devices to prune from the extracted PSY system.
- `<system_stem>_metadata.json` and `<system_stem>_validation_descriptors.json`:
  optional PSY sidecars copied with the system JSON when present.
- `simulation_settings.csv`, `markets_included.csv`, and `options.csv`: case-template
  overrides. If omitted, templates from `scripts/template/` are used.
- `queue_cost_data.csv`, investor finance defaults, and `markets_data/`: study-specific
  market and finance overrides. See
  `config/project_templates/project_spec/markets_data/README.md`.

`reference_case_dir` is optional. It can supply compatible downstream market, investor,
and finance files, but it does not replace the need for the user's own system, zones,
technologies, scenarios, options, outages, and time series.

## Generated Layout

`initialize_emis_project` writes generated project files under `output_dir`:

- `<base_dir_name>/<heterogeneity>/`: market inputs, investor inputs, and `system_config/`
- `<test_system_dir_name>/`: copied PSY system bundle
- `<test_system_dir_name>/RTS_Data/SourceData/`: derived `gen.csv`, `branch.csv`,
  `dc_branch.csv`, and `reserves.csv`
- `<test_system_dir_name>/RTS_Data/timeseries_data_files/`: copied user time series
- `case_templates/`: rendered run script, settings, market toggles, and options
- `<runs_dir_name>/`: run folders created by `new_emis_case`

`new_emis_case` only stamps a run folder from `case_templates/`. The case data folder
under `<base_dir_name>/<case_name>/` is created later by `CaseDefinition` when the run
script starts.

## Running A Case

After initialization, create a case with `new_emis_case`. Per-case settings can be
changed through overrides, for example `simulation_years=1` for a smoke run. The generated
run script copies the initialized system time series from:

```text
<test_system_dir_name>/RTS_Data/timeseries_data_files/
```

to the run result folder:

```text
<runs_dir_name>/<case_name>/Results/<run_name>/timeseries_data_files/
```

Then `create_agent_simulation` reads the generated base, test-system, market, investor,
and time series inputs.