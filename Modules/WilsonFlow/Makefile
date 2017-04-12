MKDIR=../../Make
include $(MKDIR)/MkFlags  

SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)
INCLUDE=-I ../../Include
DEPS:=$(SOURCES:%.c=%.d)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $<

%.d : %.c
	$(CC) $(MM) $(INCLUDE) -o $@ $<

default: $(OBJECTS)
$(OBJECTS): $(SOURCES)

$(DEPS): $(SOURCES)

-include $(DEPS) 

.PHONY clean:
	rm -f *.o *.d *~







