FC=gfortran

EXE1=poincare.exe
EXE2=three_bodies_legacy.exe
EXE3=three_bodies_new.exe
EXE4=three_bodies_sho.exe

all:
	$(FC) -o $(EXE1) poincare/poincare.f90 
	$(FC) -o $(EXE2) three_bodies_legacy/3_bodies_MC_legacy.f
	$(FC) -o $(EXE3) three_bodies_new/3_bodies_MC.F90 
	$(FC) -o $(EXE4) three_bodies_sho/three_bodies_sho.F90 

clean:
	rm -rf $(EXE1) $(EXE2) $(EXE3) $(EXE4) *.exe *.mod
	