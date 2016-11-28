all : simpler

#FLAGS=-xhost -O2
FC = mpif90
FLAGS = -fc=ifort 
#FLAGS=-xhost -pg -O2 (-pg is for profiling) 
#FLAGS=-g -fbacktrace #-check bounds/all (all checks everything) -fpe0 (tracks floating point errors) 

obs=decompu.o mkgridu.o calc_au.o calc_ap.o calc_aw.o mpdata_flux.o diffusion.o strainrate.o strainrate2.o fineness.o smooth.o nbrex2d2.o tracer_advec.o nbrex2d_tracer_counts.o nbrex2d_tracer_dat.o interpgen2d_geo.o nbrex2d_int.o

mods=btype.f90 bfunc2d.f90 

libs=libmg2dp.a libmgpoisson2dp.a

simpler : $(mods) $(obs) simpler.f90 $(libs) 
	$(FC) $(FLAGS) $(mods) $(obs) simpler.f90 $(libs) -o simpler 

advecdiff.o : btype.f90 bfunc2d.f90 jacobi2dfvp2.o nbrex2d.o mg2d.o advecdiff.f90 
	$(FC) $(FLAGS) -c btype.f90 bfunc2d.f90 jacobi2dfvp2.o nbrex2d.o mg2d.o advecdiff.f90

smooth.o : btype.f90 smooth.f90 
	$(FC) $(FLAGS) -c btype.f90 smooth.f90

solve_alpha.o : btype.f90 solve_alpha.f90 
	$(FC) $(FLAGS) -c btype.f90 solve_alpha.f90

fineness.o : btype.f90 fineness.f90 
	$(FC) $(FLAGS) -c btype.f90 fineness.f90  

diffusion.o : btype.f90 bfunc2d.f90 jacobi2dfvp2.o nbrex2d.o mg2d.o diffusion.f90
	$(FC) $(FLAGS) -c btype.f90 bfunc2d.f90 jacobi2dfvp2.o nbrex2d.o mg2d.o diffusion.f90 

libmg2dp.a : mg2d.o mkgrid.o decomp.o jacobi2dfvp2.o nbrex2d.o interpgen2d.o rstrct2dp.o resid2dfvp2.o 
	ar rvs libmg2dp.a mg2d.o mkgrid.o decomp.o jacobi2dfvp2.o nbrex2d.o interpgen2d.o rstrct2dp.o resid2dfvp2.o 

libmgpoisson2dp.a : mglin2d.o jacobi2dfvp2.o nbrex2d.o interpgen2d.o rstrct2dp.o resid2dfvp2.o 
	ar rvs libmgpoisson2dp.a mglin2d.o jacobi2dfvp2.o nbrex2d.o interpgen2d.o rstrct2dp.o resid2dfvp2.o 

mg2d.o : btype.f90 bfunc2d.f90 rstrct2dp.o resid2dfvp2.o interpgen2d.o nbrex2d.o jacobi2dfvp2.o mg2d.f90
	$(FC) $(FLAGS) -c btype.f90 bfunc2d.f90 rstrct2dp.o resid2dfvp2.o interpgen2d.o nbrex2d.o jacobi2dfvp2.o mg2d.f90 

mglin2d.o : btype.f90 bfunc2d.f90 rstrct2dp.o resid2dfvp.o interpgen2d.o nbrex2d.o jacobi2dfvp.o mglin2d.f90
	$(FC) $(FLAGS) -c btype.f90 bfunc2d.f90 rstrct2dp.o resid2dfvp.o interpgen2d.o nbrex2d.o jacobi2dfvp.o mglin2d.f90 

nbrex2d_tracer_dat.o : btype.f90 nbrex2d_tracer_dat.f90 
	$(FC) $(FLAGS) -c btype.f90 nbrex2d_tracer_dat.f90 

nbrex2d_tracer_counts.o : btype.f90 nbrex2d_tracer_counts.f90 
	$(FC) $(FLAGS) -c btype.f90 nbrex2d_tracer_counts.f90 

tracer_advec.o : btype.f90 tracer_advec.f90 
	$(FC) $(FLAGS) -c btype.f90 tracer_advec.f90 

strainrate2.o : btype.f90 strainrate2.f90
	$(FC) $(FLAGS) -c btype.f90 strainrate2.f90 

strainrate.o : btype.f90 strainrate.f90
	$(FC) $(FLAGS) -c btype.f90 strainrate.f90 

mpdata_flux.o : btype.f90 bfunc2d.f90 mpdata_flux.f90 
	$(FC) $(FLAGS) -c bfunc2d.f90 btype.f90 mpdata_flux.f90

#mpdata.o : btype.f90 bfunc2d.f90 mpdata.f90 
#	$(FC) $(FLAGS) -c bfunc2d.f90 btype.f90 mpdata.f90 

calc_aw.o : btype.f90 bfunc2d.f90 calc_aw.f90 
	$(FC) $(FLAGS) -c bfunc2d.f90 btype.f90 calc_aw.f90

calc_ap.o : btype.f90 bfunc2d.f90 calc_ap.f90 
	$(FC) $(FLAGS) -c bfunc2d.f90 btype.f90 calc_ap.f90

calc_au.o : btype.f90 bfunc2d.f90 calc_au.f90 
	$(FC) $(FLAGS) -c bfunc2d.f90 btype.f90 calc_au.f90

mkgridu.o : btype.f90 mkgridu.f90
	$(FC) $(FLAGS) -c btype.f90 mkgridu.f90

mkgrid.o : btype.f90 mkgrid.f90
	$(FC) $(FLAGS) -c btype.f90 mkgrid.f90

decompu.o : btype.f90 decompu.f90  
	$(FC) $(FLAGS) -c btype.f90 decompu.f90

decomp.o : btype.f90 decomp.f90  
	$(FC) $(FLAGS) -c btype.f90 decomp.f90

jacobi2dfvp2.o : btype.f90 jacobi2dfvp2.f90 
	$(FC) $(FLAGS) -c btype.f90 jacobi2dfvp2.f90 

jacobi2dfvp.o : btype.f90 jacobi2dfvp.f90 
	$(FC) $(FLAGS) -c btype.f90 jacobi2dfvp.f90 

jacobi2dp.o : btype.f90 jacobi2dp.f90 
	$(FC) $(FLAGS) -c btype.f90 jacobi2dp.f90 

nbrex2d_int.o : btype.f90 nbrex2d_int.f90 
	$(FC) $(FLAGS) -c btype.f90 nbrex2d_int.f90

nbrex2d2.o : btype.f90 nbrex2d2.f90 
	$(FC) $(FLAGS) -c btype.f90 nbrex2d2.f90

nbrex2d.o : btype.f90 nbrex2d.f90 
	$(FC) $(FLAGS) -c btype.f90 nbrex2d.f90

interpgen2d_geo.o : btype.f90 interpgen2d_geo.f90
	$(FC) $(FLAGS) -c btype.f90 interpgen2d_geo.f90

interpgen2d.o : btype.f90 interpgen2d.f90
	$(FC) $(FLAGS) -c btype.f90 interpgen2d.f90

rstrct2dp.o : btype.f90 rstrct2dp.f90 
	$(FC) $(FLAGS) -c btype.f90 rstrct2dp.f90

resid2dfvp2.o : btype.f90 resid2dfvp2.f90  
	$(FC) $(FLAGS) -c btype.f90 resid2dfvp2.f90

resid2dfvp.o : btype.f90 resid2dfvp.f90  
	$(FC) $(FLAGS) -c btype.f90 resid2dfvp.f90

clean : 
	rm *.mod *.o simpler 