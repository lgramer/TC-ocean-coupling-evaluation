#!/bin/sh

module purge
module load hpss
module list

module purge
module list
module load gnu/9.2.0 netcdf/4.7.2
module list


