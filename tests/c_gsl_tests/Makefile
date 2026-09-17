CFLAGS = -Wall -I/usr/lib -lm -lgsl -g
GRAPHFLAGS = -C

define print_cmd
	printf '  %-20s %s\n' $(1) $(2)
endef

all: out akima_plot akima_periodic_plot steffen_plot bilinear_plot bicubic_plot

out:
	@$(call print_cmd, "MKDIR", "created tests/c_gsl_tests/out/ directory")
	@mkdir out

clean:
	@$(call print_cmd, "CLEAN", "cleaned test artifacts")
	@rm -rf out

# ==============================================================

akima: out
	@$(call print_cmd, "AKIMA", "running Akima test")
	@gcc gsl_akima.c -o out/akima $(CFLAGS)
	@./out/akima > ./out/akima.dat

akima_plot: akima
	@$(call print_cmd, "AKIMA", "Saving Akima plot at tests/c_gsl_tests/out/")
	@graph -T png $(GRAPHFLAGS) < ./out/akima.dat > ./out/akima.png

# ==============================================================

akima_periodic: out
	@$(call print_cmd, "AKIMA-PERIODIC", "running Akima periodic test")
	@gcc gsl_akima_periodic.c -o out/akima_periodic $(CFLAGS)

akima_periodic_plot: akima_periodic
	@$(call print_cmd, "AKIMA-PERIODIC", "Saving Akima periodic plot at tests/c_gsl_tests/out/")
	@./out/akima_periodic > ./out/akima_periodic.dat
	@graph -T png $(GRAPHFLAGS) < ./out/akima_periodic.dat > ./out/akima_periodic.png

# ==============================================================

steffen: out
	@$(call print_cmd, "STEFFEN", "running Steffen test")
	@gcc gsl_steffen.c -o out/steffen $(CFLAGS)
	@./out/steffen > ./out/steffen.dat

steffen_plot: steffen
	@$(call print_cmd, "STEFFEN", "Saving Steffen plot at tests/c_gsl_tests/out/")
	@graph -T png $(GRAPHFLAGS) < ./out/steffen.dat > ./out/steffen.png

# ==============================================================

bilinear: out
	@$(call print_cmd, "BILINEAR", "running Bilinear test")
	@gcc gsl_bilinear.c -o out/bilinear $(CFLAGS)
	@./out/bilinear > ./out/bilinear.dat

bilinear_plot: bilinear
	@$(call print_cmd, "BILINEAR", "Saving Bilinear plot at tests/c_gsl_tests/out/")
	@$(shell ./bilinear_plot.sh)

# ==============================================================

bicubic: out
	@$(call print_cmd, "BICUBIC", "running Bicubic test")
	@gcc gsl_bicubic.c -o out/bicubic $(CFLAGS)
	@./out/bicubic > ./out/bicubic.dat

bicubic_plot: bicubic
	@$(call print_cmd, "BICUBIC", "Saving Bicubic plot at tests/c_gsl_tests/out/")
	@$(shell ./bicubic_plot.sh)
