set sty dat l;

set xlabel '(R-r)^+'
set ylabel 'U_z^+ (shifted by 5 units)'
set log x
set yrange [0:45]
set xrange [1:1200]
set key left top

set size 1.0,1.0
set term postscrip eps enhanced "Helvetica" 25 lw 2.0
set output "dns.eps"

plot \
'DNS_raw/180_Re_1.dat'  us 2:($3+0)   ti 'Re_{/Symbol t}=180'  lt 1 linecolor rgb "red", \
'DNS_raw/360_Re_1.dat'  us 2:($3+5)   ti 'Re_{/Symbol t}=360'  lt 1 linecolor rgb "green", \
'DNS_raw/550_Re_1.dat'  us 2:($3+10)  ti 'Re_{/Symbol t}=550'  lt 1 linecolor rgb "blue", \
'DNS_raw/1000_Re_1.dat' us 2:($3+15)  ti 'Re_{/Symbol t}=1000' lt 1 linecolor rgb "black"

!epstopdf dns.eps

