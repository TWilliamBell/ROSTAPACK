# Changelog

ROSTAPACK: RObust STAbility PACKage


## Version: 2.2 --- 2019-09-06

### Added

- specValSet now uses Hessenberg factorizations by default to make evaluations
  of the norm of transfer function much faster.  This results in specValSet
  being noticeably faster as well.

### Changed

- specValSet: E must be an invertible matrix
- specValSet: better tolerance for getting imaginary eigenvalues when using eig
- specValSetOptions: opts.use_permutations has been replaced by opts.solve_type
- specValSetOptions: opts.ignore_infinite has been removed

### Documentation

- Various documentation typos fixed 

### Maintenance

 - Some minor code cleanup and refactoring of private subroutines


## Version: 2.1 --- 2019-04-24

### Fixed

- Undid seemingly harmless last minute tiny bit of code cleanup which actually
  was causing numerically inaccuracy in 2nd derivatives whenever p and m weren't
  equal.  2nd derivatives should now be accurate in all cases.  In general,
  setting opts.bracket_order and opts.root_order both to 2 should be the most
  efficient configuration for specValSet.


## Version: 2.0 --- 2019-02-14

### Added

- New specValSet routine for exact computation of the spectral value set
  abscissa and radius (or the pseudospectral abscissa and radius)


## Version: 1.0.1 --- 2019-02-13

### Changed

- improved efficiency of isZero and areIdentical subroutines
- convertDSS,luSolver: renamed use_umfpack to use_permutations
- updated optionValidator to latest version

### Documentation

- corrected some typos in getStabRadBoundOptions and stateSpaceABCD
- several documentation updates for URTM subroutines
- removed last period in "included in ROSTAPACK Version 1.0" from all files,
  to avoid confusion with point releases
- moved to markdown format for CHANGELOG


## Version 1.0 --- 2018-05-18

Description: initial public release
