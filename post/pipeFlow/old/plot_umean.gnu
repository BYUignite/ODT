set sty dat l;


#------- planar means of T for given times

set log x
set yrange[0:45]
set xrange[1:1200]
set ylabel 'u^+ (shifted)';
set xlabel '(R-r)^+';
set key left top

set size 1.0,1.0
set term postscrip eps enhanced "Helvetica" 25 lw 2.0
set output "umeanRe.eps"

set arrow from 100,12 to 20,33 linecolor rgb "grey"
set label 'Re' at 100,10 textcolor "grey"

plot \
'DNS_raw/180_Re_1.dat'      us 2:3       ti 'DNS'   lt 2 linecolor rgb "blue"  , \
'ODTstats_180_2.dat'        us 1:2       ti 'ODT'   lt 1 linecolor rgb "black" , \
'DNS_raw/360_Re_1.dat'      us 2:($3+5)  ti ''      lt 2 linecolor rgb "blue"  , \
'DNS_raw/550_Re_1.dat'      us 2:($3+10) ti ''      lt 2 linecolor rgb "blue"  , \
'DNS_raw/1000_Re_1.dat'     us 2:($3+15) ti ''      lt 2 linecolor rgb "blue"  , \
'ODTstats_360_2.dat'        us 1:($2+5)  ti ''      lt 1 linecolor rgb "black" , \
'ODTstats_550_2.dat'        us 1:($2+10) ti ''      lt 1 linecolor rgb "black", \
'ODTstats_1000_2.dat'       us 1:($2+15) ti ''      lt 1 linecolor rgb "black"


! epstopdf umeanRe.eps
! rm -rf umeanRe.eps
! open -a preview umeanRe.pdf
