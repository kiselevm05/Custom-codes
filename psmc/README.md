# Scripts for PSMC analyses

## plot_psmc.py

It is an alternative to `psmc_plot.pl` included in psmc conda package. The script plots PSMC results from multiple files considering bootstrap replicates. It also allows calculating bootstrap intervals (Pseudo-
CI) and saving them to a CSV file.

usage: plot_psmc.py [-h] -f FILES [FILES ...] [-o OUTPUT] [-c [CONF_INT]] [--legend LEGEND [LEGEND ...]] [--legend-title LEGEND_TITLE]
                    [-u MUT] [-g GEN] [-x MIN_X] [-X MAX_X] [-y MIN_Y] [-Y MAX_Y] [--add-generations]

options:
  -h, --help            show this help message and exit
  
  -f FILES [FILES ...], --files FILES [FILES ...]
                        List of paths to input PSMC files (space-separated). (default: None)
                        
  -o OUTPUT, --output OUTPUT
                        Output path and filename for the plot (supports .jpg, .pdf, .png). (default: psmc_plot.jpg)
                        
  -c [CONF_INT], --conf-int [CONF_INT]
                        Calculate 95% bootstrap intervals (Pseudo-CI) and save to the specified CSV file. If no filename is provided,
                        defaults to conf-int.csv. If output path is specified, uses a directory for the output file. (default: None)
                        
  --legend LEGEND [LEGEND ...]
                        Space-separated list of legend names corresponding to the input files. Legend is displayed only if there are 2 or
                        more input files. (default: None)
                        
  --legend-title LEGEND_TITLE Legend title (default: PSMC-files)
                        
  -u MUT, --mut MUT     Mutation rate per nucleotide per generation. (default: 2.5e-08)
  
  -g GEN, --gen GEN     Generation time in years. (default: 25)
  
  -x MIN_X, --min-x MIN_X
                        Minimum X-axis value in years. (default: 10000.0)
                        
  -X MAX_X, --max-x MAX_X
                        Maximum X-axis value in years. (default: 10000000.0)
  -y MIN_Y, --min-y MIN_Y
                        Minimum Y-axis value (x10^4 of individuals). (default: 0)
                        
  -Y MAX_Y, --max-y MAX_Y
                        Maximum Y-axis value (x10^4 of individuals). (default: None)
                        
  --add-generations     Display the number of generations on the X-axis in parentheses after the year value. (default: False)`

  ## compare_psmc.py

Compare Ne values of two PSMC assemblies at a specific time point using bootstrap replicates and the Mann-Whitney U test.

  usage: compare_psmc.py [-h] -f1 FILE1 -f2 FILE2 -t TIME [-u MUT] [-g GEN]

options:
  -h, --help            show this help message and exit
  
  -f1 FILE1, --file1 FILE1  First .psmc file (default: None)
                        
  -f2 FILE2, --file2 FILE2  Second .psmc file (default: None)
  
  -t TIME, --time TIME  Target time in years before present, integer (default: None)
  
  -u MUT, --mut MUT     Mutation rate per site per generation (-u in psmc_plot) (default: 2.5e-08)
  
  -g GEN, --gen GEN     Generation time (-g in psmc_plot) (default: 25)
