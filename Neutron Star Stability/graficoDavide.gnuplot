



# Set terminal to PNG or PDF for high-quality output
# set terminal pngcairo size 800,600 enhanced font 'Helvetica,12'
#set output 'graficoerrorerkeulero.png'

set tmargin at screen 0.9
# Set title and labels
set grid
set lmargin at screen 0.15
set rmargin at screen 0.95

set xlabel '' font ",14"
set ylabel '' font ",14"
set title "" font ",16"

# Set plot style and colors
# Add more styles if needed
# Set grid
set grid

set logscale y
set logscale x


# Plot data
plot "./data/filegraficodavide.dat" using 3:1 w lp ps 1 pt 5 notitle, \
"./data/filegraficodavide.dat" using 3:2 w lp ps 1 pt 5 notitle,
# Add more data sets as needed
# Add legend
set key top left

# Add any other customizations
