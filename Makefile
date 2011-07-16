TASK = star_stat_db


include /proj/sot/ska/include/Makefile.FLIGHT
SHARE = update_star_stats.py agasc.py
DATA = ok_missing.csv task_schedule.cfg


install: 
ifdef BIN
	mkdir -p $(INSTALL_BIN)
	rsync --times --cvs-exclude $(BIN) $(INSTALL_BIN)/
endif
ifdef SHARE
	mkdir -p $(INSTALL_SHARE)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)
endif
ifdef DATA
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)
endif
ifdef LIB
	mkdir -p $(INSTALL_PERLLIB)/$(PERLGEN)
	rsync --times --cvs-exclude $(LIB) $(INSTALL_PERLLIB)/$(PERLGEN)/
endif
ifdef WWW
	mkdir -p $(SKA)/www/ASPECT/acq_stats/
	rsync --times --cvs-exclude $(WWW) $(SKA)/www/ASPECT/acq_stats/
endif
