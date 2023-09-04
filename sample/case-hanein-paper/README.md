# Simulation of cases in reference paper

This folder replicates some cases evaluated by [Hanein (2017)](https://doi.org/10.1080/17436753.2017.1303261). The supplementary data from their paper is reconsolidated under [data/](data/) directory, with both experimental conditions and results in a single file for each of their sources (44 cases from Tscheng and 9 cases from Barr). In both cases thermocouple temperatures are given in Kelvin and coordinate of 0 meters represents the solid feed end of the kiln.

**Note:** because the model does not currently accept zero-thickness coating, an arbitrary value of 3 mm was added (and removed from refractory), the value also being used to correct internal radius (`diameter = 2 * radius + 6` to account coating over diameter). It is possible to verify that this procedure ensures `model.R_cv` provides the right value of internal radius reported in the paper.

## Detailed data from Barr

- Solid loading is set at 12%.

- Kiln rotation rate is set at 1.5 rpm.

- Gas temperature at 10 cm off the kiln wall labeled `Tg_off_wall`.

- Gas temperature at 2.5 cm off the kiln bed labeled `Tg_off_bed`.

- Solid bed and wall temperatures are labeled `Ts` and `Tw`, respectively.

- The numbers following column headers denote the thermocouple numbers.

    1. The gas thermocouple locations are: (1) 0.11, (2) 0.89, (3) 2.15, (4) 2.51, (5) 2.85, (6) 3.20, (7) 3.91, (8) 4.44, and (9) 4.95.

    1. The solid bed thermocouple locations are: (1) 0.11, (2) 0.87, (3) 1.44, (4) 2.14, (5) 2.50, (6) 2.85, (7) 3.20, (8) 3.91, (9) 4.48, (10) 4.95, (11) 5.25, and (12) 5.50 meters.

    1. The wall thermocouple locations are: (1) 1.33, (2) 2.32, (3) 2.67, (4) 3.03, (5) 3.38, (6) 3.83, (7) 4.40, and (8) 4.78.

## Detailed data from Tscheng

- Solid particle diameters are 0.73 mm.  

- Temperatures are labelled `T-M-#` where `M` denotes the material in question (`g` for gas, `s` for solid bed, and `w` for wall) and `#` denotes the thermocouple number.

- The locations of the gas and bed thermocouples are: (1) 0.21, (2) 0.72, (3) 1.25, (4) 1.78, and (5) 2.32 meters along the kiln.

- The locations of the wall thermocouples are: (1) 0.31, (2) 0.91, (3) 1.52, and (4) 2.13 meters along the kiln.
