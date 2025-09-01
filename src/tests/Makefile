CFLAGS = -Wall -I/usr/lib -lm -lgsl -g

all: out akima_plot akima_periodic_plot steffen_plot

out:
	mkdir out -v

clean:
	rm -rf out

# ==============================================================

akima: out
	gcc gsl_akima.c -o out/akima $(CFLAGS)
	@./out/akima > ./out/akima.dat

akima_plot: akima
	@graph -T png < ./out/akima.dat > ./out/akima.png 
	$(info Plot saved at src/gsl_tests/out/akima.png)

# ==============================================================

akima_periodic: out
	gcc gsl_akima_periodic.c -o out/akima_periodic $(CFLAGS)

akima_periodic_plot: akima_periodic
	@./out/akima_periodic > ./out/akima_periodic.dat
	@graph -T png < ./out/akima_periodic.dat > ./out/akima_periodic.png 
	$(info Plot saved at src/gsl_tests/out/akima_periodic.png)

# ==============================================================

steffen: out
	gcc gsl_steffen.c -o out/steffen $(CFLAGS)
	@./out/steffen > ./out/steffen.dat

steffen_plot: steffen
	@graph -T png < ./out/steffen.dat > ./out/steffen.png 
	$(info Plot saved at src/gsl_tests/out/steffen.png)
