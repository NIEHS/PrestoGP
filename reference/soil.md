# Soil data from North Bavaria, Germany

Soil physical and chemical data collected on a field in the
Weissenstaedter Becken, Germany

## Usage

``` r
data(soil)
```

## Format

A data frame with the following columns:

- x.coord::

  x coordinate (cm)

- y.coord::

  y coordinate (cm)

- nr::

  Sample number. The order of the sample number corresponds to the order
  that the samples were collected.

- moisture::

  Moisture content (Kg/Kg \* 100%)

- NO3.N::

  Nitrate nitrogen (mg/Kg)

- Total.N::

  Total nitrogen (mg/Kg)

- NH4.N::

  Ammonium nitrogen (mg/Kg)

- DOC::

  Dissolved organic carbon (mg/Kg)

- N20N::

  Nitrous oxide (mg/Kg dried substance)

## Source

The data were collected by Wolfgang Falk, Soil Physics Group, University
of Bayreuth, Germany. This data set was previously published in the
"RandomFields" R package.

## Details

For technical reasons some of the data were obtained as differences of
two measurements (which are not available anymore). Therefore, some of
the data have negative values.

## References

- Falk, W. "Kleinskalige raeumliche Variabilitaet von Lachgas und
  bodenchemischen Parameters (Small scale spatial variability of nitrous
  oxide and pedo-chemical parameters)" (2000). Master thesis, University
  of Bayreuth, Germany.

## Author

Martin Schlather, School of Business Informatics and Mathematics,
University of Mannheim
