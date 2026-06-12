# Exoplanet Observing Window Planner

A Python command-line tool for calculating observable transit and eclipse windows for exoplanet observations.

This script was developed to support high-resolution exoplanet spectroscopy planning, where observing windows must satisfy multiple practical constraints including target visibility, airmass, twilight, Moon separation, exposure time, and instrument overheads.

The tool queries planet parameters, evaluates observability from a specified observing site, optionally generates diagnostic plots, and outputs timing windows in a format that can be used for observation planning.

## Features

* Calculates observable transit windows for exoplanet targets
* Calculates pre-eclipse and post-eclipse observing windows
* Applies practical observing constraints including:

  * nautical twilight
  * maximum airmass
  * Moon separation
  * exposure time
  * instrument overhead
  * required observing-window length
* Queries planet parameters using the NASA Exoplanet Archive
* Uses `astroplan` and `astropy` for target visibility and timing calculations
* Uses `batman` to model transit and eclipse phases
* Generates optional diagnostic plots showing orbital phase and airmass
* Outputs timing windows in an observation-planning-friendly text format

## Research Context

High-resolution exoplanet spectroscopy often requires precise scheduling around transits, eclipses, and orbital phases. Because these observations are sensitive to target visibility, atmospheric conditions, and instrument constraints, planning requires combining astrophysical ephemerides with practical observability limits.

This tool automates that process by identifying candidate observing windows that satisfy both scientific and operational constraints.

## Installation

Clone the repository:

```bash
git clone https://github.com/emilydeibert/exoplanet-observing-window-planner.git
cd exoplanet-observing-window-planner
```

Install the required packages:

```bash
pip install -r requirements.txt
```

## Requirements

This tool uses:

```text
astroplan
astroquery
astropy
matplotlib
numpy
batman-package
```

## Usage

The script can be run from the command line. For example, to calculate observable transit windows:

```bash
python observing_window_planner.py --planet "WASP-121 b" --type transit --exp_time 300 --plot
```

To calculate eclipse observing windows:

```bash
python observing_window_planner.py --planet "WASP-121 b" --type eclipse --exp_time 300 --length 5 --airmass 2 --plot
```

## Command-Line Arguments

| Argument              | Description                                                               |
| --------------------- | ------------------------------------------------------------------------- |
| `-p`, `--planet`      | Planet name formatted for the NASA Exoplanet Archive, e.g. `"WASP-121 b"` |
| `-t`, `--type`        | Type of observing window to calculate: `transit` or `eclipse`             |
| `-plot`, `--plot`     | Generate optional diagnostic plots                                        |
| `-e`, `--exp_time`    | Exposure time in seconds                                                  |
| `-length`, `--length` | Minimum eclipse observing-window length in hours; default is 5            |
| `-a`, `--airmass`     | Maximum allowed airmass; default is 2                                     |

## Example Output

The tool writes timing-window files with the extension:

```text
.tw
```

These files contain observing-window start times and durations in a format intended to be useful for observation planning.

Example output format:

```text
YYYY-MM-DD HH:MM:SS HH:MM
```

where the final column gives the observing-window duration.

## Example Workflow

```bash
python observing_window_planner.py --planet "WASP-121 b" --type transit --exp_time 300 --airmass 2 --plot
```

This command will:

1. Query planet parameters for WASP-121 b.
2. Identify transit windows during the configured observing semester.
3. Apply observability constraints.
4. Add available timing flexibility where possible.
5. Save the resulting timing windows to a `.tw` file.
6. Generate diagnostic plots if `--plot` is enabled.

## Notes and Limitations

This is research software developed for exoplanet observation planning. Some parameters, including observing site, semester dates, instrument overheads, and special-case planet parameters, are currently configured directly in the script.

Future improvements could include:

* Moving observatory, semester, and instrument settings into a configuration file
* Adding more general support for multiple observatories
* Improving handling of missing or uncertain archive parameters
* Adding automated tests for timing and formatting functions
* Refactoring special-case planet parameters into a separate data file
