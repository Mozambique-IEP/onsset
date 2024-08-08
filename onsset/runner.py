# Defines the modules

import logging
import os
import time

import geojson
import pandas as pd
from onsset import (SET_ELEC_ORDER, SET_LCOE_GRID, SET_MIN_GRID_DIST, SET_GRID_PENALTY,
                    SET_MV_CONNECT_DIST, SET_WINDVEL, SET_WINDCF, SettlementProcessor, Technology)

try:
    from onsset.specs import (SPE_COUNTRY, SPE_ELEC, SPE_ELEC_MODELLED,
                              SPE_ELEC_RURAL, SPE_ELEC_URBAN, SPE_END_YEAR,
                              SPE_GRID_CAPACITY_INVESTMENT, SPE_GRID_LOSSES,
                              SPE_MAX_GRID_EXTENSION_DIST,
                              SPE_NUM_PEOPLE_PER_HH_RURAL,
                              SPE_NUM_PEOPLE_PER_HH_URBAN, SPE_POP, SPE_POP_FUTURE,
                              SPE_START_YEAR, SPE_URBAN, SPE_URBAN_FUTURE,
                              SPE_URBAN_MODELLED)
except ImportError:
    from specs import (SPE_COUNTRY, SPE_ELEC, SPE_ELEC_MODELLED,
                       SPE_ELEC_RURAL, SPE_ELEC_URBAN, SPE_END_YEAR,
                       SPE_GRID_CAPACITY_INVESTMENT, SPE_GRID_LOSSES,
                       SPE_MAX_GRID_EXTENSION_DIST,
                       SPE_NUM_PEOPLE_PER_HH_RURAL,
                       SPE_NUM_PEOPLE_PER_HH_URBAN, SPE_POP, SPE_POP_FUTURE,
                       SPE_START_YEAR, SPE_URBAN, SPE_URBAN_FUTURE,
                       SPE_URBAN_MODELLED)
from openpyxl import load_workbook

logging.basicConfig(format='%(asctime)s\t\t%(message)s', level=logging.DEBUG)


def calibration(specs_path, csv_path, specs_path_calib, calibrated_csv_path):
    """

    Arguments
    ---------
    specs_path
    csv_path
    specs_path_calib
    calibrated_csv_path
    """
    specs_data = pd.read_excel(specs_path, sheet_name='SpecsData')
    settlements_in_csv = csv_path
    settlements_out_csv = calibrated_csv_path

    onsseter = SettlementProcessor(settlements_in_csv)

    num_people_per_hh_rural = float(specs_data.iloc[0][SPE_NUM_PEOPLE_PER_HH_RURAL])
    num_people_per_hh_urban = float(specs_data.iloc[0][SPE_NUM_PEOPLE_PER_HH_URBAN])

    # RUN_PARAM: these are the annual household electricity targets
    tier_1 = 38.7  # 38.7 refers to kWh/household/year. It is the mean value between Tier 1 and Tier 2
    tier_2 = 219
    tier_3 = 803
    tier_4 = 2117
    tier_5 = 3000

    onsseter.prepare_wtf_tier_columns(tier_1, tier_2, tier_3, tier_4, tier_5)
    onsseter.condition_df()
    onsseter.df[SET_GRID_PENALTY] = onsseter.grid_penalties(onsseter.df)

    onsseter.df[SET_WINDCF] = onsseter.calc_wind_cfs(onsseter.df[SET_WINDVEL])

    pop_actual = specs_data.loc[0, SPE_POP]
    urban_current = specs_data.loc[0, SPE_URBAN]
    start_year = int(specs_data.loc[0, SPE_START_YEAR])
    elec_actual = specs_data.loc[0, SPE_ELEC]
    elec_actual_urban = specs_data.loc[0, SPE_ELEC_URBAN]
    elec_actual_rural = specs_data.loc[0, SPE_ELEC_RURAL]

    pop_modelled, urban_modelled = onsseter.calibrate_current_pop_and_urban(pop_actual, urban_current)

    specs_data.loc[0, SPE_URBAN_MODELLED] = urban_modelled

    elec_calibration_results = onsseter.calibrate_grid_elec_current(elec_actual, elec_actual_urban, elec_actual_rural,
                                                                    start_year, buffer=False)

    mg_pop_calib = onsseter.mg_elec_current(start_year)

    specs_data.loc[0, SPE_ELEC_MODELLED] = elec_calibration_results[0]
    specs_data.loc[0, 'rural_elec_ratio_modelled'] = elec_calibration_results[1]
    specs_data.loc[0, 'urban_elec_ratio_modelled'] = elec_calibration_results[2]
    specs_data['grid_distance_used'] = elec_calibration_results[3]
    specs_data['ntl_limit'] = elec_calibration_results[4]
    specs_data['pop_limit'] = elec_calibration_results[5]
    specs_data['Buffer_used'] = elec_calibration_results[6]
    specs_data['buffer_distance'] = elec_calibration_results[7]
    specs_data['mg_pop_electrified'] = mg_pop_calib

    book = load_workbook(specs_path)
    with pd.ExcelWriter(specs_path_calib, engine='openpyxl') as writer:
        writer.workbook = book
        # RUN_PARAM: Here the calibrated "specs" data are copied to a new tab called "SpecsDataCalib".
        # This is what will later on be used to feed the model
        specs_data.to_excel(writer, sheet_name='SpecsDataCalib', index=False)

    logging.info('Calibration finished. Results are transferred to the csv file')
    onsseter.df.to_csv(settlements_out_csv, index=False)


def scenario(specs_path, calibrated_csv_path, results_folder, summary_folder, pv_path):
    """

    Arguments
    ---------
    specs_path : str
    calibrated_csv_path : str
    results_folder : str
    summary_folder : str

    """

    scenario_info = pd.read_excel(specs_path, sheet_name='ScenarioInfo')
    scenarios = scenario_info['Scenario']
    scenario_parameters = pd.read_excel(specs_path, sheet_name='ScenarioParameters')
    specs_data = pd.read_excel(specs_path, sheet_name='SpecsDataCalib', index_col=0)
    print(specs_data.iloc[0][SPE_COUNTRY])

    for scenario in scenarios:
        print('Scenario: ' + str(scenario + 1))

        yearsofanalysis = specs_data.index.tolist()
        base_year = specs_data.iloc[0][SPE_START_YEAR]
        end_year = yearsofanalysis[-1]
        start_years = [base_year] + yearsofanalysis

        time_steps = {}
        for year in range(len(yearsofanalysis)):
            time_steps[yearsofanalysis[year]] = yearsofanalysis[year] - start_years[year]

        onsseter = SettlementProcessor(calibrated_csv_path)


        x_coordinates, y_coordinates = onsseter.start_extension_points(r'C:\Users\andre\OneDrive\Dokument\GitHub\SEforALL-onsset\test_data\MV_lines_guess.geojson')
        onsseter.add_xy_3395()

        country_id = specs_data.iloc[0]['CountryCode']

        # ToDo project pop for each year based on previous year
        pop_future = specs_data.iloc[0][SPE_POP_FUTURE]
        urban_future = specs_data.iloc[0][SPE_URBAN_FUTURE]

        # ToDo make more flexible to read all from one sheet ???
        tier_index = scenario_info.iloc[scenario]['Target_electricity_consumption_level']
        grid_index = scenario_info.iloc[scenario]['Grid_electricity_generation_cost']
        pv_index = scenario_info.iloc[scenario]['PV_cost_adjust']
        productive_index = scenario_info.iloc[scenario]['Productive_uses_demand']
        prio_index = scenario_info.iloc[scenario]['Prioritization_algorithm']

        rural_tier = scenario_parameters.iloc[tier_index]['RuralTargetTier']
        urban_tier = scenario_parameters.iloc[tier_index]['UrbanTargetTier']
        grid_price = scenario_parameters.iloc[grid_index]['GridGenerationCost']
        pv_capital_cost_adjust = scenario_parameters.iloc[pv_index]['PV_Cost_adjust']
        productive_demand = scenario_parameters.iloc[productive_index]['ProductiveDemand']
        prioritization = scenario_parameters.iloc[prio_index]['PrioritizationAlgorithm']
        auto_intensification = scenario_parameters.iloc[prio_index]['AutoIntensificationKM']
        max_auto_intensification_cost = scenario_parameters.iloc[prio_index]['MaxIntensificationCost']  # Max household connection cost for forced grid intensification

        settlements_out_csv = os.path.join(results_folder,
                                           '{}-1-{}_{}_{}_{}.csv'.format(country_id, tier_index, grid_index, pv_index,
                                                                         prio_index))
        summary_csv = os.path.join(summary_folder,
                                   '{}-1-{}_{}_{}_{}_summary.csv'.format(country_id, tier_index, grid_index, pv_index,
                                                                         prio_index))

        elements = ["1.Population", "2.New_Connections", "3.Capacity", "4.Investment", "5.AnnualEmissions"]

        techs = ["Grid", "SA_PV", "MG_Diesel", "MG_PVHybrid", "MG_Wind", "MG_Hydro"]
        tech_codes = [1, 3, 4, 5, 6, 7]

        sumtechs = []
        for element in elements:
            for tech in techs:
                sumtechs.append(element + "_" + tech)
        total_rows = len(sumtechs)
        df_summary = pd.DataFrame(columns=yearsofanalysis)
        for row in range(0, total_rows):
            df_summary.loc[sumtechs[row]] = "Nan"

        onsseter.current_mv_line_dist()

        onsseter.project_pop_and_urban(pop_future, urban_future, base_year, yearsofanalysis)

        discount_rate = 0.08  # RUN_PARAM

        carbon_cost = scenario_parameters['CarbonTax'][grid_index]  # ToDo does not need to be linked to grid index?
        grid_emission_factor = scenario_parameters['GridEmissionFactor'][grid_index]
        grid_capacity_investment = scenario_parameters['GridCapacityInvestmentCost'][grid_index]
        grid_re_share = scenario_parameters['GridRenShare'][grid_index]


        # Carbon cost represents the cost in USD/tonCO2eq, which is converted and added to the diesel price
        diesel_price = float(scenario_parameters.iloc[0]['DieselPrice'] + (carbon_cost / 1000000) * 256.9131097 * 9.9445485)

        for year in yearsofanalysis:

            time_step = time_steps[year]
            start_year = year - time_step

            # RUN_PARAM: Fill in general and technology specific parameters (e.g. discount rate, losses etc.)
            Technology.set_default_values(base_year=start_year,
                                          start_year=start_year,
                                          end_year=end_year,
                                          discount_rate=discount_rate)

            grid_calc = Technology(om_of_td_lines=0.02,
                                   distribution_losses=float(specs_data.iloc[0][SPE_GRID_LOSSES]),
                                   connection_cost_per_hh=125,
                                   base_to_peak_load_ratio=0.8,
                                   capacity_factor=1,
                                   tech_life=30,
                                   grid_capacity_investment=grid_capacity_investment,
                                   grid_penalty_ratio=1,
                                   grid_price=grid_price)

            mg_hydro_calc = Technology(om_of_td_lines=0.02,
                                       distribution_losses=0.05,
                                       connection_cost_per_hh=100,
                                       base_to_peak_load_ratio=0.85,
                                       capacity_factor=0.5,
                                       tech_life=30,
                                       capital_cost={float("inf"): 3000},
                                       om_costs=0.03,
                                       mini_grid=True)

            mg_wind_calc = Technology(om_of_td_lines=0.02,
                                      distribution_losses=0.05,
                                      connection_cost_per_hh=100,
                                      base_to_peak_load_ratio=0.85,
                                      capital_cost={float("inf"): 3750},
                                      om_costs=0.02,
                                      tech_life=20,
                                      mini_grid=True)

            # mg_pv_calc = Technology(om_of_td_lines=0.02,
            #                         distribution_losses=0.05,
            #                         connection_cost_per_hh=100,
            #                         base_to_peak_load_ratio=0.85,
            #                         tech_life=20,
            #                         om_costs=0.015,
            #                         capital_cost={float("inf"): 2950 * pv_capital_cost_adjust},
            #                         mini_grid=True)

            sa_pv_calc = Technology(base_to_peak_load_ratio=0.9,
                                    tech_life=15,
                                    om_costs=0.02,
                                    capital_cost={float("inf"): 6950 * pv_capital_cost_adjust,
                                                  1: 4470 * pv_capital_cost_adjust,
                                                  0.100: 6380 * pv_capital_cost_adjust,
                                                  0.050: 8780 * pv_capital_cost_adjust,
                                                  0.020: 9620 * pv_capital_cost_adjust
                                                  },
                                    standalone=True)

            mg_interconnection = True  # True if mini-grids are allowed to be integrated into the grid, else False

            mg_hybrid_lookup_table = True
            mg_pv_hybrid_params = {
                'min_mg_size_ppl': 500,  # minimum number of people in settlement for mini-grids to be considered as an option
                'diesel_cost': 261,  # diesel generator capital cost, USD/kW rated power
                'discount_rate': discount_rate,
                'n_chg': 0.93,  # charge efficiency of battery
                'n_dis': 1,  # discharge efficiency of battery
                'battery_cost': 314,  # battery capital cost, USD/kWh of storage capacity
                'pv_cost': 660,  # PV panel capital cost, USD/kW peak power
                'charge_controller': 142,  # PV charge controller cost, USD/kW peak power, set to 0 if already included in pv_cost
                'pv_inverter': 80,  # PV inverter cost, USD/kW peak power, set to 0 if already included in pv_cost
                'pv_life': 25,  # PV panel expected lifetime, years
                'diesel_life': 10,  # diesel generator expected lifetime, years
                'pv_om': 0.015,  # annual OM cost of PV panels
                'diesel_om': 0.1,  # annual OM cost of diesel generator
                'battery_inverter_cost': 539,
                'battery_inverter_life': 20,
                'dod_max': 0.8,  # maximum depth of discharge of battery
                'inv_eff': 0.93,  # inverter_efficiency
                'lpsp_max': 0.02,  # maximum loss of load allowed over the year, in share of kWh
                'diesel_limit': 0.5,  # Max annual share of mini-grid generation from diesel gen-set
                'full_life_cycles': 2500  # Equivalent full life-cycles of battery until replacement
            }

            mg_diesel_calc = Technology(om_of_td_lines=0.02,
                                        distribution_losses=0.05,
                                        connection_cost_per_hh=100,
                                        base_to_peak_load_ratio=0.85,
                                        capacity_factor=0.7,
                                        tech_life=15,
                                        om_costs=0.1,
                                        capital_cost={float("inf"): 721},
                                        mini_grid=True)

            sa_diesel_calc = Technology(base_to_peak_load_ratio=0.9,
                                        capacity_factor=0.5,
                                        tech_life=10,
                                        om_costs=0.1,
                                        capital_cost={float("inf"): 938},
                                        standalone=True)

            sa_diesel_cost = {'diesel_price': diesel_price,
                              'efficiency': 0.28,
                              'diesel_truck_consumption': 14,
                              'diesel_truck_volume': 300}

            mg_diesel_cost = {'diesel_price': diesel_price,
                              'efficiency': 0.33,
                              'diesel_truck_consumption': 33.7,
                              'diesel_truck_volume': 15000}

            eleclimit = specs_data.loc[year]['ElecTarget']
            num_people_per_hh_rural = float(specs_data.loc[year][SPE_NUM_PEOPLE_PER_HH_RURAL])
            num_people_per_hh_urban = float(specs_data.loc[year][SPE_NUM_PEOPLE_PER_HH_URBAN])
            max_grid_extension_dist = float(specs_data.loc[year][SPE_MAX_GRID_EXTENSION_DIST])
            annual_grid_cap_gen_limit = specs_data.loc[year, 'NewGridGenerationCapacityAnnualLimitMW'] * 1000 * time_step
            annual_new_grid_connections_limit = specs_data.loc[year]['GridConnectionsLimitThousands'] * 1000 * time_step

            onsseter.calculate_demand(year, num_people_per_hh_rural, num_people_per_hh_urban, time_step,
                                      urban_tier, rural_tier)

            onsseter.diesel_cost_columns(sa_diesel_cost, mg_diesel_cost, year)

            if mg_hybrid_lookup_table:
                hybrid_lcoe, hybrid_capacity, hybrid_investment = \
                    onsseter.pv_hybrids_lcoe_lookuptable(year, time_step, end_year,
                                                         mg_pv_hybrid_params, pv_path=pv_path)
            else:
                hybrid_lcoe, hybrid_capacity, hybrid_investment = \
                    onsseter.pv_hybrids_lcoe(year, time_step, end_year,
                                             mg_pv_hybrid_params, pv_folder_path=pv_path)

            mg_pv_hybrid_calc = Technology(om_of_td_lines=0.02,
                                           distribution_losses=0.05,
                                           connection_cost_per_hh=100,
                                           capacity_factor=0.5,
                                           base_to_peak_load_ratio=0.85,
                                           tech_life=20,
                                           mini_grid=True,
                                           hybrid_fuel=hybrid_lcoe,
                                           hybrid_investment=hybrid_investment,
                                           hybrid_capacity=hybrid_capacity,
                                           hybrid=True)

            sa_diesel_investment, sa_diesel_capacity, sa_pv_investment, sa_pv_capacity, mg_diesel_investment, \
            mg_diesel_capacity, mg_pv_hybrid_investment, mg_pv_hybrid_capacity, mg_wind_investment, mg_wind_capacity, \
            mg_hydro_investment, mg_hydro_capacity = onsseter.calculate_off_grid_lcoes(mg_hydro_calc, mg_wind_calc,
                                                                                       sa_pv_calc,
                                                                                       mg_diesel_calc, sa_diesel_calc,
                                                                                       mg_pv_hybrid_calc,
                                                                                       year, end_year, time_step,
                                                                                       techs, tech_codes)

            grid_investment, grid_capacity, grid_cap_gen_limit, grid_connect_limit = \
                onsseter.pre_electrification(grid_price, year, time_step, end_year, grid_calc,
                                             annual_grid_cap_gen_limit, annual_new_grid_connections_limit)

            onsseter.max_extension_dist(year, time_step, end_year, start_year, grid_calc, max_auto_intensification_cost)

            onsseter.pre_selection(eleclimit, year, time_step, prioritization, auto_intensification)

            new_lines_geojson = {}
            onsseter.df[SET_LCOE_GRID + "{}".format(year)], onsseter.df[SET_MIN_GRID_DIST + "{}".format(year)], \
                grid_investment, grid_capacity, x_coordinates, y_coordinates, new_lines_geojson[year] = \
                onsseter.elec_extension_numba(grid_calc,
                                              max_grid_extension_dist,
                                              year,
                                              start_year,
                                              end_year,
                                              time_step,
                                              grid_cap_gen_limit,
                                              grid_connect_limit,
                                              grid_investment,
                                              grid_capacity,
                                              x_coordinates,
                                              y_coordinates,
                                              auto_intensification=auto_intensification,
                                              prioritization=prioritization,
                                              threshold=max_auto_intensification_cost,
                                              )

            onsseter.results_columns(techs, tech_codes, year, time_step, prioritization, auto_intensification,
                                     mg_interconnection)

            onsseter.calculate_investments_and_capacity(sa_diesel_investment, sa_diesel_capacity, sa_pv_investment,
                                                        sa_pv_capacity, mg_diesel_investment, mg_diesel_capacity,
                                                        mg_pv_hybrid_investment, mg_pv_hybrid_capacity, mg_wind_investment,
                                                        mg_wind_capacity, mg_hydro_investment, mg_hydro_capacity,
                                                        grid_investment, grid_capacity, year)

            if year == yearsofanalysis[-1]:
                final_step = True
            else:
                final_step = False

            onsseter.check_grid_limitations(annual_new_grid_connections_limit, annual_grid_cap_gen_limit, year, time_step, final_step)

            onsseter.apply_limitations(eleclimit, year, time_step, prioritization, auto_intensification)

            onsseter.calculate_emission(grid_factor=grid_emission_factor, year=year,
                                        time_step=time_step, start_year=start_year)

            onsseter.calc_summaries(df_summary, sumtechs, tech_codes, year, base_year)

            # Save to a GeoJSON file
            with open(os.path.join(results_folder, 'new_mv_lines_{}_{}.geojson'.format(scenario, year)), 'w') as f: # ToDo
                geojson.dump(new_lines_geojson[year], f)

        for i in range(len(onsseter.df.columns)):
            if onsseter.df.iloc[:, i].dtype == 'float64':
                onsseter.df.iloc[:, i] = pd.to_numeric(onsseter.df.iloc[:, i], downcast='float')
            elif onsseter.df.iloc[:, i].dtype == 'int64':
                onsseter.df.iloc[:, i] = pd.to_numeric(onsseter.df.iloc[:, i], downcast='signed')

        df_summary.to_csv(summary_csv, index=sumtechs)
        onsseter.df.to_csv(settlements_out_csv, index=False)

        logging.info('Finished')
