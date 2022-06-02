set output 'pf1000.eps'
set terminal postscript eps color enhanced font "Helvetica,20"
set xlabel "t [{/Symbol ns}]"
set ylabel "I"
plot "test.out" u 1:2 with points title "I"
