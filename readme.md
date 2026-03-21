## Table of contents

1. [Description of the repository]
2. [Requirements]
3. [Repository constitution]
4. [References]

# Description of the repository

This repository is devoted to the determination of interfacial energy in binary silicate systems with a miscibility gap. The interfacial energy is determined from the model developed by Kaptay [1]. Using an alternative method to compute the interfacial energy according to Cahn and Hilliard [2], the interfacial thickness is computed.

The computations are based on the thermodynamic data describing the thermodynamic equilibrium of the separeted phases. The Gibbs free energy is determined using the excess potential of Redlich-Kister potential [3]. Eleven systems are considered for which coefficients of the Redlich-Kister potential are gathered in a file, dbexcesspotential.csv. The experimental data of the critical and monotectic conditions are in critcondbinarysystem.csv and monotecticbinarysystem.csv files.

The details of the development of the models can be found in the submitted article:

F. Pigeonneau, I. Martin and W. Blanc (2026). "Interfacial energy and thickness in phase separated silicate systems" J. Am. Ceram. Soc.

## Requirements

Python modules involved in this repository are:

-   `numpy`
-   `pandas`
-   `matplotlib`
-   `scipy`
-    'molarmass'

## Repository constitution

Provided python files are:

-   'binarysystem.py'
-   'interfaciallayer.py'
-   'phasediagram.py'
-   'surfacetensiontwocomponent.py'

Data coming from experimental data and for Redlich-Kister potential are in the following csv files:

-   'critcondbinarysystem.csv'
-   'dbexcesspotential.csv'
-   'dboxides.csv'
-   'monotecticbinarysystem.csv'

Experimental data described the binodal bounds for all systems are gathered in the directory 'BinodalData'. Molar fractions as a function of the temperature are given in files 'xvsT-Systemname.csv' for which 'Systemname' correspond to the name of system studied in this repository. Eleven systems are provided for this work.

## References

[1] G. Kaptay (2012). "On the interfacial energy of coherent interfaces". Acta Mater., 60:6804-6813.
[2] J. W. Cahn & J. E. Hilliard (1958). "Free Energy of a Nonuniform System". I. Interfacial Free Energy. J. Chem. Phys., 28:258-267.
[3] O. Redlich & A. T. Kister (1948). "Algebraic representation of thermodynamic properties and the classification of solutions". Ind. Eng. Chem., 40:345-348.
