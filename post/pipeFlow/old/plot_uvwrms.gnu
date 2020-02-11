set sty dat l;


#------- planar means of T for given times

set yrange[0:3.5]
set xrange[-1000:1000]
set ylabel 'u_{RMS}^+'
set xlabel '(R-r)^+';
set key left top

set size 1.0,1.0
set term postscrip eps enhanced "Helvetica" 25 lw 2.0
set output "uvwRms1000.eps"

set arrow nohead from 0,0 to 0,3.5 linecolor rgb "grey"
set arrow        from  400,1.7 to  700,1.7 linecolor rgb "grey"
set arrow        from -400,1.7 to -700,1.7 linecolor rgb "grey"
set label 'ODT' at  400,2 textcolor rgb "grey"
set label 'DNS' at -600,2 textcolor rgb "grey"

plot \
'ODTstats_1000_2.dat'           us 1:5 ti ''           lt 1 linecolor rgb "black", \
'ODTstats_1000_2.dat'           us 1:6 ti ''           lt 3 linecolor rgb "blue",  \
'ODTstats_1000_2.dat'           us 1:7 ti ''           lt 5 linecolor rgb "red",  \
'DNS_raw/DNS_1000_reversed.dat' us 2:7 ti 'u_{z, RMS}' lt 1  linecolor rgb "black",  \
'DNS_raw/DNS_1000_reversed.dat' us 2:5 ti 'u_{r, RMS}' lt 3  linecolor rgb "blue",  \
'DNS_raw/DNS_1000_reversed.dat' us 2:6 ti 'u_{t, RMS}' lt 5  linecolor rgb "red"

! epstopdf uvwRms1000.eps
! rm -rf uvwRms1000.eps
! open -a preview uvwRms1000.pdf
