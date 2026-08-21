## What this module does

This repository contains the source code of the Open Source Spatial Electrification Tool (OnSSET), adapted to the case of Mozambique. This includes codes to calibrate the base year information, as well as to run geospatial electrification scenarios to identify the least-cost technology (grid, mini-grid and Solar Home System) for each settlement in the country, as well as associated capacity and investment requirements.

## Installation

The extraction module (as well as all supporting scripts in this repo) have been developed in Python 3. We recommend installing [Anaconda's free distribution](https://www.anaconda.com/distribution/) as suited for your operating system. In order to be able to run the tool you have to install all necessary packages contained in "moz_onsset_env.yml". To do this, simply open Anaconda prompt and browse to the code directory on your computer and run:

```
conda env create --name moz_onsset_env --file moz_onsset_env.yml
```

## How to run

1. Download input data to data/inputs (see paths in **Input data** section below)
2. Install the required packages (see installation instructions using Anaconda above)
3. Activate the environment in Anaconda prompt (*conda activate moz_onsset_env*)
4. Launch Jupyter Notebook (*jupyter notebook* in Anaconda Prompt)
5. There are two notebooks to be run:
  * a) **1_Calibration** is run to prepare the start year information (identification of electrified locations, calibration of population to match official statistics, etc.)
  * b) **2_Scenarios** is run to calculate the least-cost technology options for the target year. 
  * If you already have a calibrated file, you can proceed directly to the **2_Scenarios** notebook, otherwise first run **1_Calibration**.
  The third notebook in the repository is called **detailed_tech_parameters**. Here one can specify additional input parameters for the scenario runs, but the most common and important input parameters are typically found in the scenarios notebook.
6. Choose the notebook and follow the instructions inside
7. Outputs will be written to data/outputs
8. Check that data in the output csv file match with expected values as described in **Output data** below
9. If you run **1_Calibration** to generate the calibrated input file, you must move or copy the file from the output folder to the correct input folder in order to use it for the scenario runs.


## Input data

### Calibration

| Dataset | Target location | Description | SDI Location | Alternative location |
|---------|-----------------|-------------|--------------|----|
| InputFile.csv | data/inputs/extracted_csv | Input file with extracted GIS-data for every settlement in Mozambique | https://piem.mireme.gov.mz/en/admin/datasets/referencedataset/6 | - |

### Scenario run(s)

| Dataset | Target location | Description | SDI Location | Alternative location |
|---------|-----------------|-------------|--------------|----|
| Calibrated.csv | data/inputs/calibrated_csv | Calibrated input file with GIS-data and calibrated data for every settlement in Mozambique | https://piem.mireme.gov.mz/en/admin/datasets/referencedataset/7 | - |
| MV_existing_and_proposed.geojson | data/inputs/MV | GIS-file of existing MV lines as well as any lines already committed to be built | https://piem.mireme.gov.mz/en/admin/datasets/vectordataset/32/ | - |
| mz-2-pv.csv | data/inputs/pv | CSV-file with hourly PV data for Mozambique | https://piem.mireme.gov.mz/en/admin/datasets/referencedataset/8 | [EnergyData.Info](https://energydata.info/dataset/mozambique-global-electrification-platform-gep) |
| mz-2-wind.csv | data/inputs/wind | CSV-file with hourly wind data for Mozambique | https://piem.mireme.gov.mz/en/admin/datasets/referencedataset/9 | [EnergyData.Info](https://energydata.info/dataset/mozambique-global-electrification-platform-gep) |

## Output data

### Calibration

The output data is a file with the extracted GIS data for every settlement, as well as the newly added calibrated information, saved as a csv (*OnSSET_InputFile_Calibrated.csv*) in outputs/calibration. 
The csv-file contains all of the columns described in the Output data section of the [OnSSET-GIS-Extraction repo](https://github.com/Mozambique-IEP/OnSSET-GIS-Extraction/tree/master#output-data), as well as the following **new** columns:

| Column | Type | Unit | Description |
|---------|-----------------|-------------|--------------|
| PopStartYear |	Float |	people |	Population in the settlement in the start year of the analysis |
| ElecPopCalib |	Float |	people |	Calibrated version of ElecPop in the start of the analysis |
| NumPeoplePerHH |	Float |	people |	Number of people in household; value for urban or rural respective to how the settlement is characterized |
| Pop20*XX* |	Float |	people |	Projected population at the year |
| ElecStart |	Integer |	0-1 |	0 if settlement not electrified today; 1 if electrified today; Retrieved from calibration algorithm |
| FinalElecCode20*XX* |	Integer |	0-99 |	Code defining type of technology providing electricity if electrified in base year |
| ElecPop20*XX* |	Float |	people |	Population in settlement assumed to have access to electricity in the specific year |
| HHs20*XX* | Float | households | Households in the settlement in the start year of the analysis |

### Scenario(s)

The first output data is a file with extracted GIS data for every settlement, saved as a csv (*Scenario_name*_Results.csv* - Note that ***Scenario_name*** can be replaced with any name by the user) in outputs/scenarios:
This file contains all of the columns included in the *OnSSET_InputFile_Calibrated.csv* file, as well as:

| Column | Type | Unit | Description |
|---------|-----------------|-------------|--------------|
| WindCF |	Float |	%	| Estimated capacity factor for wind technologies; certain technical parameters are taken into consideration |
| Pop20*XX* |	Float |	people |	Projected population at the year |
| NewConnections20*XX* |	Float |	Households |	Households that need to be electrified by the year |
| EnergyPerSettlement20*XX* |	Float |	kWh |	Estimated electricity demand target in the settlement based on new connections |
| MG_Hydro20*XX* |	Float |	USD/kWh |	LCOE estimated for mini-grid hydro in the year |
| MG_PVHybrid20*XX* |	Float |	USD/kWh |	LCOE estimated for mini-grid pv/diesel/battery hybrid in the year |
| MG_Wind20*XX* |	Float |	USD/kWh |	LCOE estimated for mini-grid wind/diesel/battery hybrid in the year |
| SA_PV20*XX* |	Float |	USD/kWh |	LCOE estimated for stand-alone PV (SHS) in the year |
| Minimum_Tech_Off_grid20*XX* |	String |	- |	The off-grid technology that can meet the demand at the lowest LCOE in the year |
| Minimum_LCOE_Off_grid20*XX* |	Float |	USD/kWh |	LCOE estimated for the minimum off-grid technology in the year |
| OffGridInvestmentCost20*XX* |	Float |	USD |	Investment cost for the minimum off-grid technology in the year |
| Grid20*XX* |	Float |	USD/kWh |	LCOE estimated for the Grid in the year |
| NewGridExtensionDist20*XX* |	Float |	km |	Indicating the distance of the grid, if settlement gets grid electrified in the process |
| MinimumOverall20*XX* |	String |	tech abbreviation |	Abbreviation defining type of off-grid technology providing electricity in the step year |
| MinimumOverallLCOE20*XX* |	Float |	USD/kWh |	LCOE of the least cost off grid option selected |
| MinimumOverallCode20*XX* |	Integer |	1 - 7 code |	Code defining type of off-grid technology providing electricity in the step year |
| InvestmentCost20*XX* |	Float |	USD |	Total investment if electrification is achieved |
| InvestmentPerConnection20*XX* |	Float |	USD/capita |	Estimated investment per connection if electrification is achieved |
| ElecStatusIn20*XX* |	Integer |	0-1 |	Final electrification status in the step year (after running prioritization algorithm) |
| NewCapacity20*XX* |	Float |	kW |	Additional capacity for the least cost technology identified in the step year |
| TotalEnergyPerCell |	Float |	kWh/year |	Total electricity demand (all sectors considered) in the settlement by the end year of the analysis |
| Tier |	Int |	1-5 |	Tier classification of the consumption per capita in the settlement |
| Technology20*XX* |	String |	- |	Name of the technology in the settlement in the year, or ‘unelectrified’ |
| AnnualEmissions20*XX* |	Float |	kg CO2 eq. / year |	Estimated new emissions induced by the selected technology for new connections in that year |

The second output data is a file with summary results for the entire country, saved as a csv (*Scenario_name*_Summaries.csv* - Note that ***Scenario_name*** can be replaced with any name by the user) in outputs/scenarios.
The following summary information is provided for each technology type as well as combined (in the Total row):

| Column | Type | Unit | Description |
|---------|-----------------|-------------|--------------|
| Population20*XX* |	Float |	people |	Population served by each technology (or Non-electrified) by that year |
| NewConnections20*XX* |	Float |	households |	Households connected by each technology in that timestep |
| Capacity20*XX* |	Float |	MW |	New generation capacity added by each technology in that timestep |
| Investment20*XX* |	Float |	Million USD |	Investment required by each technology in that timestep |
| Emissions20*XX* |	Float |	tonne CO2 eq. |	Annual new emissions from each technology in that timestep |
| PopulationTotal |	Float |	people |	Population served by each technology (or Non-electrified) by the end of the analysis |
| NewConnectionsTotal |	Float |	households |	Households connected by each technology in both time-steps combined |
| CapacityTotal |	Float |	MW |	New generation capacity added by each technology in both time-steps combined |
| InvestmentTotal |	Float |	Million USD |	Investment required by each technology in both time-steps combined |
| EmissionsTotal |	Float |	tonne CO2 eq. |	Annual new emissions from each technology in both time-steps combined |

The third output data is a csv file with key input parameters used to define the scenario, saved as a csv (*Scenario_name*_Variables.csv* - Note that ***Scenario_name*** can be replaced with any name by the user) in outputs/scenarios. 

The fourth output data is a GIS file of new MV lines to extend the grid in the scenario, which can be used for visualization purposes. It is saved as a geojson (*Scenario_name*_new_mv_lines.geojson* - Note that ***Scenario_name*** can be replaced with any name by the user) in outputs/scenarios. 

## License

This module is made available under the **MIT** license.
See the LICENSE file in this repository for the full text.

## Contact

Questions/maintainer:
- Inocencio Gujamo - MIREME-UIPCE (inocencio.gujamo@gmail.com)
- Imaculada Dos Santos - MIREME-UIPCE (imaculadamz@gmail.com)

Developed by: 
- Andreas Sahlberg - SEforALL (andreas.sahlberg@seforall.org)
- Alexandros Korkovelos - SEforALL (alexandros.korkovelos@seforall.org)
- Babak Khavari - SEforALL (babak.khavari@seforall.org)
- Julian Cantor - SEforALL (julian.cantor@seforall.org)

