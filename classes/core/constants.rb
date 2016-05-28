################################################################################
#
# Constants
#
# Author: Jan Rezac
# Date created: 2008-11-17
# License: Cuby license
# Description: Unit conversion factors and constants
# Status: Works, Documented
#
################################################################################

#: Units used in Cuby

#: Fundamental units:
#* length 	- angstrom
#* energy 	- kcal/mol
#* time   	- picosecond
#* temperature 	- kelvin
#* charge	- elementary charge (charge of a proton)

#: Derived units:
#* velocity  	- A/ps = 100 m/s :-)
#* mass   	- kcal/mol.ps^2/A^2 = 418.4 g/mol
#* force  	- kcal/mol/A

################################################################################
# Elementary constants from NIST
################################################################################

#===============================================================================
# Charge - conversion factors
#===============================================================================

#: from http://physics.nist.gov/cuu/Constants/index.html
ELEMENTARY2COULOMB = 1.602176487e-19
COULOMB2ELEMENTARY = 1.0/1.602176487e-19

#===============================================================================
# Universal constants in SI units
#===============================================================================
#: from http://physics.nist.gov/cuu/Constants/index.html
SPEED_OF_LIGHT_SI 	= 299792458.0		# m s-1
PLANCK_CONSTANT_SI 	= 6.62606896e-34	# J s
ELEMENTARY_CHARGE_SI	= 1.602176487e-19	# C
AVOGADRO_CONSTANT_SI	= 6.02214179e23		# mol-1
BOLTZMANN_CONSTANT_SI	= 1.3806504e-23		# J K-1
GAS_CONSTANT_SI		= 8.314472		# J mol-1 K-1
ELECTRON_REST_MASS_SI	= 9.10938215e-31 	# kg

# Shortcuts
AVOGADRO = AVOGADRO_CONSTANT_SI

#===============================================================================
# Atomic units in SI
#===============================================================================

# This section is to be removed
HARTREE_SI		= 4.35974394e-18 	# J
BOHR_SI			= 0.52917720859e-10 	# m


################################################################################
# Unit conversion
################################################################################

# Length
ANGSTROM2BOHR = (1.0/0.5291772083)
BOHR2ANGSTROM = 0.5291772083

NANOMETER2ANGSTROM = 10.0
ANGSTROM2NANOMETER = 0.1

ANGSTROM2M = 1.0e-10
M2ANGSTROM = 1.0e10

# Energy
HARTREE2KCAL = 627.5095
KCAL2HARTREE = (1.0/627.5095)

CAL2J = 4.184
J2CAL = 1.0/CAL2J

KCAL2KJ =4.184
KJ2KCAL = 1.0/KCAL2KJ

KCAL2J = CAL2J * 1000.0
J2KCAL = 1.0/CAL2J/1000.0

J2EV = 1.0/ELEMENTARY_CHARGE_SI
EV2J = 1.0/J2EV

# Warning: used as is, it produces EV/mol
EV2KCAL = 23.060538
KCAL2EV = 1.0/23.060538

EV2HARTREE = EV2KCAL * KCAL2HARTREE
HARTREE2EV = 1.0 / EV2HARTREE

# Mass
GMOL2UNIT = (1.0/418.4)
UNIT2GMOL = (418.4)

MASSUNIT2KG = UNIT2GMOL / 1000.0 / AVOGADRO

# Time
PS2S = 1.0e-12

# Frequency
HZ2CM = 1.0 / 100.0 / SPEED_OF_LIGHT_SI # Hz -> cm-1, 1/100/c
CM2HZ = 1.0 / HZ2CM

CM2KCALMOL = 100.0 * SPEED_OF_LIGHT_SI * PLANCK_CONSTANT_SI * J2KCAL * AVOGADRO
KCALMOL2CM = 1.0 / CM2KCALMOL

# Pressure
ATM2PA = 101325.0

# Dipole
# Dimension Length * Charge
DEBYE2AU = 0.393430307
AU2DEBYE = 1.0 / DEBYE2AU
DEBYE2CUBY = DEBYE2AU * BOHR2ANGSTROM
CUBY2DEBYE = 1.0 / DEBYE2CUBY


################################################################################
# Constants in our units
################################################################################

BOLTZMANN = 0.0019872 # kcal.mol-1.K-1 
PLANCK_CONSTANT = PLANCK_CONSTANT_SI * J2KCAL * 1.0e12 # kcal ps


################################################################################
# Atomic units
################################################################################

# Constants that are equall to 1.0:
# * elementary charge
# * electron rest mass
# * Bohr radius
# * abs. value of electronic enegy of H atom in ground state
# * Planck constant
# * Constant in Coulomb's law: 1.0/(4*pi*epsilon0)

# Conversion from SI
#: from http://physics.nist.gov/cuu/Constants/index.html
AU2COULOMB	= ELEMENTARY_CHARGE_SI	# Charge
AU2KG		= ELECTRON_REST_MASS_SI	# Mass
AU2M		= 0.52917720859e-10	# Length
AU2J		= 4.35974394e-18	# Energy
# Derived units
AU2S		= 2.418884326505e-17	# Time
AU2M_S		= 2.1876912541e6	# Velocity
AU2NEWTON	= 8.23872206e-8		# Force
AU2K		= 3.1577464e5		# Temperature, from http://en.wikipedia.org/wiki/Atomic_units

# Reverse
COULOMB2AU	= 1.0 / AU2COULOMB
KG2AU		= 1.0 / AU2KG
M2AU		= 1.0 / AU2M
J2AU		= 1.0 / AU2J
S2AU		= 1.0 / AU2S
M_S2AU		= 1.0 / AU2M_S
NEWTON2AU	= 1.0 / AU2NEWTON
K2AU		= 1.0 / AU2K

#===============================================================================
# Conversion from CUBY units
# Charge - the same unit
MASS2AU		= UNIT2GMOL / 1000.0 * KG2AU / AVOGADRO
ANGSTROM2AU	= ANGSTROM2BOHR
KCAL2AU		= KCAL2HARTREE	# KCAL means kcal/mol
PS2AU		= 1.0e-12 * S2AU
CUBYVELO2AU	= ANGSTROM2AU / PS2AU
CUBYFORCE2AU	= MASS2AU * ANGSTROM2AU / PS2AU**2

#Reverse
AU2MASS		= 1.0 / MASS2AU
AU2ANGSTROM	= 1.0 / ANGSTROM2AU
AU2KCAL		= 1.0 / KCAL2AU
AU2PS		= 1.0 / PS2AU
AU2CUBYVELO	= 1.0 / CUBYVELO2AU
AU2CUBYFORCE	= 1.0 / CUBYFORCE2AU

#===============================================================================
# Derived constants in AU
BOLTZMANN_AU = BOLTZMANN_CONSTANT_SI * J2AU / K2AU
