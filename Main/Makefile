MKDIR=../Make
include $(MKDIR)/MkFlags
PROF= #-pg

MOD=../Modules/
EXE = TWI_Z_exact TWI_Z_WFJac main_const_phys main_coupling main_TWI_RC main_DWI_RC main_TDWI_RC main_new_coupling
SRC= TWI_Z_exact.c TWI_Z_WFJac.c main_const_phys.c main_coupling.c main_TWI_RC.c main_DWI_RC.c main_TDWI_RC.c main_new_coupling.c
OBJ = $(SRC:.c=.o)
DEPS=$(SRC:.c=.d)
INC=-I ../Include

DIRS=$(sort $(dir $(wildcard $(MOD)*/)))
MODOBJ=$(wildcard $(MOD)*/*.o)


all:dirs $(EXE)

.PHONY: dirs $(DIRS)

dirs: $(DIRS)

$(DIRS):
	$(MAKE) -C $@ 


$(OBJ): $(SRC)

$(DEPS): $(SRC)

-include $(DEPS)

$(EXE): $(OBJ)

.PHONY cleandir:
	-for d in $(DIRS); do (cd $(MOD)/$$d; $(MAKE) clean ); done

clean:cleandir
	rm -f $(EXE) *.o *~ *.d

%.o : %.c
	$(CC) $(CFLAGS) $(INC) -o $@ $<

%.d : %.c
	$(CC) $(MM) $(INC) -o $@ $<

% : %.o 
	$(CC) $(MODOBJ) $< $(INC) -lm -o $@ $(PROF) 