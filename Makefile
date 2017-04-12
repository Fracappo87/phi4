#Some conventions#
CC=gcc
CFLAG=  -Wall -O3  -msse -msse2 -pg -g
CFLAGS= -c  -Wall  -O3  -msse -msse2 -pg -g
LIST= geometry.o input_output.o lattice_config.o update.o observables.o ranlxd.o Wilson_Flow.o Jacobian.o testRK.o testJACOB.o


all: main_const_phys main_crt_point main_testRK main_testJACOB main_testMASS main_TWI_jac_ex main_corr main_TWI_sq

input_output.o: input_output.c Makefile
	 $(CC) $(CFLAGS) input_output.c

geometry.o: geometry.c Makefile
	$(CC) $(CFLAGS) geometry.c 

lattice_config.o: lattice_config.c Makefile
	$(CC) $(CFLAGS) lattice_config.c

update.o: update.c Makefile
	$(CC) $(CFLAGS) update.c

observables.o: observables.c Makefile
	$(CC) $(CFLAGS) observables.c

ranlxd.o: ranlxd.c Makefile
	$(CC) $(CFLAGS) ranlxd.c

Wilson_Flow.o: Wilson_Flow.c Makefile
	$(CC) $(CFLAGS) Wilson_Flow.c

Jacobian.o: Jacobian.c Makefile
	$(CC) $(CFLAGS) Jacobian.c

testRK.o: testRK.c Makefile
	$(CC) $(CFLAGS) testRK.c

testJACOB.o: testRK.c Makefile
	$(CC) $(CFLAGS) testJACOB.c

main_const_phys: main_const_phys.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_const_phys   main_const_phys.c $(LIST) -lm   

main_crt_point: main_crt_point.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_crt_point   main_crt_point.c $(LIST) -lm   

main_testRK: main_testRK.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_testRK   main_testRK.c $(LIST) -lm  

main_testJACOB: main_testJACOB.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_testJACOB   main_testJACOB.c $(LIST) -lm 

main_testMASS: main_testMASS.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_testMASS   main_testMASS.c $(LIST) -lm # -lefence

main_TWI_jac_ex:main_TWI_jac_ex.c  $(LIST) Makefile
	$(CC) $(CFLAG)  -o main_TWI_jac_ex   main_TWI_jac_ex.c $(LIST) -lm # -lefence

main_corr: main_corr.c $(LIST) Makefile
	$(CC) $(CFLAG) -o main_corr   main_corr.c $(LIST) -lm 

main_TWI_sq : main_TWI_sq.c $(LIST) Makefile
	$(CC) $(CFLAG) -o main_TWI_sq   main_TWI_sq.c $(LIST) -lm 

clean:
	rm -rf $(LIST) *~
