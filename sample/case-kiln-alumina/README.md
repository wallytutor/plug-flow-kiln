# Alumina Kiln

## Hypothesis

- Alumina properties for alpha variant. Data extracted from JANAF as reported [here](https://webbook.nist.gov/cbook.cgi?ID=C1344281&Type=JANAF&Table=on). Since there is a single range reported up to 2327 K, but the program requires low/high ranges, the same Shomate polynomials are used in both. Commented-out code is provided in [simulate.m](./simulate.m) for verifying this values against the source.

- Since alumina grain size is below the minimum value tested by Li (2005), a saturation air film thickness of 0.25 (dimensionless) is chosen.

- Residual O2 at fumes outlet is 3%v (used for computing discharge secondary air).
