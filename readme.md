# Table of contents

1. [Description of the repository]
2. [Requirements]
3. [Repository constitution]
4. [References]

## Description of the repository

This repository is devoted to the determination of interfacial energy in binary silicate and borate systems with a miscibility gap. Excess potentials are determined according to the Redlich and Kister model up to the fourth order [1]. The interfacial energy is determined from the model developed by Kaptay [2].

The computations are based on the thermodynamic data describing the thermodynamic equilibrium of the separeted phases. The Gibbs free energy is determined using the excess potential of Redlich-Kister potential [3]. Fourteen systems are considered for which coefficients of the Redlich-Kister potential are gathered in a file, dbexcesspotential.csv. The experimental data of the critical and monotectic conditions are in critcondbinarysystem.csv and monotecticbinarysystem.csv files.

The details of the development of the models can be found in the reference [4].

## Requirements

Python modules involved in this repository are:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `molarmass`

## Repository constitution

Provided python files are:

-   `binarysystem.py`: BinarySystem class developed to determine thermodynamic functions and their derivatives of binary system. Methods to find equilibrium state, critical condition, monotectic points and interfacial energy are provided.
-   `RedlichKisteroptimisation.py`: python script using the BinarySystem class to determine the coefficients of the Redlich-Kister function by a minimization of the difference between experimental and model predictions.
-   `surfacetensiontwocomponent.py`: python script using the BinarySystem class to determine the interfacial energy for binary system.

Data coming from experimental data and for Redlich-Kister potential are in the following csv files:

-   `critcondbinarysystem.csv`: file gathering the critical conditions from experimental result for fourteen binary silicate systems.
-   `dbexcesspotential.csv`: file gathering the coefficients of the Redlich-Kister potential at the fourth order coming from literature or determined for this contribution. 
-   `dboxides.csv`: file gathering the useful properties of oxides (density, fusion temperature and enthalpy).
-   `monotecticbinarysystem.csv`: file gathering monotectic conditions of bynary silicate systems. 

Experimental data described the binodal bounds for all systems are gathered in the directory 'BinodalData'. Molar fractions as a function of the temperature are given in files 'Systemname.csv' for which 'Systemname' correspond to the name of system studied in this repository. Fourteen systems are provided for this work.

## References

[1] O. Redlich & A. T. Kister (1948). "Algebraic representation of thermodynamic properties and the classification of solutions". Ind. Eng. Chem., 40:345-348.

[2] G. Kaptay (2012). "On the interfacial energy of coherent interfaces". Acta Mater., 60:6804-6813.

[3] J. W. Cahn & J. E. Hilliard (1958). "Free Energy of a Nonuniform System". I. Interfacial Free Energy. J. Chem. Phys., 28:258-267.

[4] F. Pigeonneau and W. Blanc (2026). "Interfacial energy in phase separated borate and silicate systems". J. Am. Ceram. Soc., 10.1111/JACE.71109.
