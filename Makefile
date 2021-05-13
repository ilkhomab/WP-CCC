#  e = new.o scat-in.o zvode_lapack.o bess_wp.o bessel.o intrpl.o filon.o qp_H_structure.o phase.o wignerd.o cwf8.o HContinuumGautchi.o coul90.o main.o
  e = wignerd.o functions_dev.o gaulag.o intrpl.o qp_H_structure.o cwf0.o HContinuumGautchi.o coul90.o iqpack.o main.o

  a = -fast -mp -r4 -acc -ta=tesla:cc60
  a = -g -mp -r4 -acc -ta=tesla:managed -traceback -Minfo=accel
  a = -g -C -Ktrap=ovf -mp -acc -ta=tesla:cc60 -traceback -Minfo=accel
  f = -g -C -mp -acc -ta=tesla:cuda10.1 -traceback -Minfo=accel
  a = -fast -mp -acc -Minfo=all -ta=tesla:cuda10.1 -mcmodel=medium -llapack -lblas -Mcudalib=cusolver -Mcuda=cuda10.1 -g77libs
#  a = -fast -mp -acc -Minfo=all -ta=tesla:cc35 -mcmodel=medium -llapack -lblas -Mcudalib=cusolver -Mcuda=cc35 -g77libs
#  a = -g -mp -r4 -acc -ta=tesla:cc30 -traceback -Minfo=accel




#b = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -reentrancy threaded
#b = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -reentrancy threaded
#  b = -Wl,--start-group -Wl,--end-group -openmp# -reentrancy threaded
#b = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64
#b = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -reentrancy threaded
b = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -mp 


#d = -L/group/d35/ilkhom/shiv -lmap-sampler-pmpi -lmap-sampler -Wl,--eh-frame-hdr -Wl,-rpath=/group/d35/ilkhom/shiv 

#  LD =  mpif90   #ifort
#  LD =  ifort 
  LD =  mpif90 

  
main_dev: $e 
	$(LD) $a $e -o main_dev

zvode_lapack.o: zvode_lapack.f
	$(LD) $a -c zvode_lapack.f

scat.o: scat.f90
	$(LD) $a -c scat.f90

gaulag.o: gaulag.f90
	gfortran -O3  -c gaulag.f90

iqpack.o: iqpack.f
	gfortran -O3  -c iqpack.f

scat-in.o: scat-in.f90
	$(LD) $a -c scat-in.f90

main.o: main.f90
	$(LD) $a -c  main.f90

bessel.o: bessel.f
	$(LD) $a -c bessel.f

bess_wp.o: bess_wp.f
	$(LD) $a -c bess_wp.f

new.o: new.f90
	$(LD) $a -c new.f90

dev.o: dev.f90
	$(LD) $a -c dev.f90

dev2.o: dev2.f90
	$(LD) $a -c dev2.f90

functions_dev.o: functions_dev.f90
	$(LD) $a -c functions_dev.f90

functions.o: functions.f90
	$(LD) $a -c functions.f90

wigd.o: wigd.f90
	$(LD) $a -c wigd.f90

wignerd.o: wignerd.f90
	$(LD) $a -c wignerd.f90

phase.o: phase.f90
	$(LD) $a -c phase.f90

filon.o: filon.f90
	$(LD) $a -c filon.f90

cwf0.o: cwf0.f90
	$(LD) $a -c cwf0.f90

cwf2.o: cwf2.f90
	$(LD) $a -c cwf2.f90

cwf3.o: cwf3.f90
	$(LD) $a -c cwf3.f90

cwf4.o: cwf4.f90
	$(LD) $a -c cwf4.f90

cwf5.o: cwf5.f90
	$(LD) $a -c cwf5.f90

cwf6.o: cwf6.f90
	$(LD) $a -c cwf6.f90

cwf7.o: cwf7.f90
	$(LD) $a -c cwf7.f90

cwf8.o: cwf8.f90
	$(LD) $a -c cwf8.f90

cwf9.o: cwf9.f90
	$(LD) $a -c cwf9.f90

cwf10.o: cwf10.f90
	$(LD) $a -c cwf10.f90
intrpl.o: intrpl.f
	$(LD) $a -c intrpl.f

HContinuumGautchi.o: HContinuumGautchi.F 
	$(LD) $a -c HContinuumGautchi.F

coul90.o: coul90.f 
	$(LD) $a -c coul90.f

qp_H_structure.o: qp_H_structure.f90
	$(LD) $a -c qp_H_structure.f90
