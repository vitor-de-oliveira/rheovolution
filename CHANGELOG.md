# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [1.0.1] - 2026-10-15

### Added
- Data skip parameter as user input.
- ODE integrator scheme as user input from the options available at GSL's ``gsl_odeiv2.h``.

### Changed
- User input for longitude of the ascending node from OMEGA(deg) to Omega(deg).
- Automatic timestep to a tenth of the shortest timescale instead of the exact shortest timescale.

### Fixed
- Memory leak and heap overflow in the input file parser.
- Calculation of the angle between spin axis and figure axis and between body's lowest moment of inertia and its relative position vector.

## [1.0.0] - 2024-12-12

### Added
- Initial stable release.