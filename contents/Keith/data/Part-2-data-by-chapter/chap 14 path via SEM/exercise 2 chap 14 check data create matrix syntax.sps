* Encoding: windows-1252.

DATASET ACTIVATE DataSet3.
RECODE race (MISSING=SYSMIS) (4=0) (ELSE=1) INTO Minority.
VARIABLE LABELS  Minority 'Minority v White'.
Value labels Minority 0 'White' 1 'Minority'. 
EXECUTE.

CORRELATIONS
  /VARIABLES=Minority,byses,bytests,f1s36a2,ffugrad
   /PRINT=SIG
  /STATISTICS DESCRIPTIVES
  /matrix out(*)
  /MISSING=LISTWISE .

Frequencies
  /VARIABLES=Minority,f1s36a2
   /STATISTICS=DEFAULT.

DESCRIPTIVES Minority,byses,bytests,f1s36a2,ffugrad  
   /STATISTICS=MIN Max mean stddev
   /MISSING=LISTWISE.


