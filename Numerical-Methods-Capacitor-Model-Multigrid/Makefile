# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
L9 = /usr/lib/libblas.dylib
L10 = /usr/lib/liblapack.dylib
L11 = /usr/local/gfortran/lib/libgfortran.dylib
L12 = ../SRC
LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) $(L9) $(L10) -lm

COPT = -O3 -Wall

agmgdir=../SRC

libc= $(L11)

list_c= $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o


default: main

clean: 
	rm *.o 
	rm main

main: main.c tableauPotentiel.o affPotentiel.o prob.o time.o umfpack.o residu.o rho.o agmgSolveur.o GS.o restriction.o prolongation.o multigrille.o nombreElem.o resVect.o residuEvolution.o conjGrad.o
	cd $(agmgdir);make dseq
	cc $(COPT) $^ -o $@ $(LIB) $(list_c) $(libc)

umfpack.o: umfpack.c
	cc $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include


%.o: %.c
	cc $(COPT) -c $< -o $@


