AM_CPPFLAGS = `pkg-config --cflags libxml-2.0` $(extra_inc)  
AM_CFLAGS = @CFLAGS@ -g -Wall 

#bin_PROGRAMS =  rse rse2 rse3 prod kerneldensity tomo
bin_PROGRAMS =  rse rse2 rse3 prod kerneldensity


#############################
rse_SOURCES = lsqr.c lsqr.h \
	           lsqr_wrapper.c lsqr_wrapper.h \
	           catch_sig.c rse.c 

rse_LDFLAGS = `pkg-config --libs libxml-2.0`
rse_LDADD = -lm -lsparse $(extra_ld)

#############################
rse2_SOURCES = lsqr.c lsqr.h \
	           lsqr_wrapper.c lsqr_wrapper.h \
	           catch_sig.c rse2.c 

rse2_LDFLAGS = `pkg-config --libs libxml-2.0`
rse2_LDADD = -lm -lsparse $(extra_ld)

################################
rse3_SOURCES = lsqr.c lsqr.h \
	           lsqr_wrapper.c lsqr_wrapper.h \
	           catch_sig.c rse3.c 

rse3_LDFLAGS = `pkg-config --libs libxml-2.0`
rse3_LDADD = -lm -lsparse $(extra_ld)


##################################
kerneldensity_SOURCES = lsqr.c lsqr.h \
                   lsqr_wrapper.c lsqr_wrapper.h \
                   catch_sig.c kernel_density.c

kerneldensity_LDFLAGS = `pkg-config --libs libxml-2.0`
kerneldensity_LDADD = -lm -lsparse $(extra_ld)


##################################
prod_SOURCES = lsqr.c lsqr.h \
	           lsqr_wrapper.c lsqr_wrapper.h \
	           catch_sig.c prod.c 

prod_LDFLAGS = `pkg-config --libs libxml-2.0`
prod_LDADD = -lm -lsparse $(extra_ld)


##################################
tomo_SOURCES = export.c export.h \
	      lsqr.c lsqr.h \
	      options.h \
	      extern.h \
	      parse_options.c parse_options.h \
	      lsqr_wrapper.c lsqr_wrapper.h \
	      irm.c irm.h \
	      relative_tt.c relative_tt.h \
	      regul.c regul.h \
	      catch_sig.c catch_sig.h \
		  tomo.c

tomo_LDFLAGS = `pkg-config --libs libxml-2.0`
tomo_LDADD = -lm -lpopt -lxml2 -lmesh -lraydescartes -lsparse $(extra_ld)

#############################

