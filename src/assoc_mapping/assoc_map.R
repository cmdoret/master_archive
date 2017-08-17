# In this script I use a modified version of the basic test for allelic association
# to test for association with homozygosity in males. The test should be performed 
# only within a category of family (inferred from diploid male production).
# Cyril Matthey-Doret
# 18.08.2017

# I will use the following contingency table:
# Observed (O): 
# Males Females Total
#--------------------
# Mhom | Fhom | Thom
# Mhet | Fhet | Thet
# Mtot | Ftot | Ttot

# Expected (E):
#       Males             Females       Total
#-----------------------------------------------
# (Thom*Mtot)/Ttot | (Thom*Ftot)/Ttot | Thom
# (Thet*Mtot)/Ttot | (Thet*Ftot)/Ttot | Thet
#        Mtot      |       Mtot       | Ttot
# And measure association's chi-square as follows:
# X2 = sum{(Ox-Ex)/Ex} where x is Mhom, Mhet, Fhom, Fhet