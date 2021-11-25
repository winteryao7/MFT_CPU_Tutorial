FLAGS=-O
MY_BIN=.

SUBS = sacio.o
PROGRAM = gsact lsac2 sac_e saclst

all: $(PROGRAM) sliding_wfcc_fix_v5 XmengStackShift

$(PROGRAM): %:%.o $(SUBS)
	$(LINK.c) -o $(MY_BIN)/$@ $@.o $(SUBS) -lm

sliding_wfcc_fix_v5: sliding_wfcc_fix_v5.o
	$(LINK.c) -o $(MY_BIN)/$@ $@.o -lm

XmengStackShift: XmengStackShift.o
	$(LINK.c) -o $(MY_BIN)/$@ $@.o -lm
	
backup:
	tar cvf WFCC_src.tar *
	gzip WFCC_src.tar

clean:
	rm *.o 
