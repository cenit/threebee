#!/gnuplot
FILE9="Initial_conditions.dat" 
FILE19="Evoluted.dat" 
FILE20="Evoluted_with_noise.dat" 
FILE_SHO4="fort.4"
FILE_SHO7="fort.7"
FILE_A="Area.dat"
FILE_SG="Sole_Giove.dat" 
set terminal png enhanced 
set output 'Plot.png'
set xlabel "x'" 
set ylabel "v_x" 
set xrange   [-1:1]  #[0.4:0.9]   # [-2 :2]          #  [0.4:0.9]   
set yrange   [-1:1] #[-0.4:0.4]  # [-2 :2]          #  [-0.4:0.4]
# set arrow from  -1.2,0 to 1.2,0   lw 1   lc rgb "black"  nohead
# set arrow from  0,-1.2 to 0,1.2   lw 1   lc rgb "black"  nohead 
#set logscale y
plot FILE_SHO7 u  1:2  w   points   pt  5   lc rgb "cyan"    ps   2   t  'x  v_x'   ,\
     FILE9     u  2:4  w   points   pt  5   lc rgb "blue"    ps  0.5  t  'x  v_x'   ,\
     FILE19    u  2:4  w   points   pt  5   lc rgb "yellow"  ps  0.1  t  'x  v_x'   ,\
     FILE20    u  2:4  w   points   pt  5   lc rgb "purple"  ps  0.2  t  'x  v_x'   ,\
     FILE_SHO4 u  2:4  w   points   pt  5   lc rgb "green"   ps  0.2  t  'x  v_x'   ,\
     FILE_A    u  2:3  w   lines    lt  1   lc rgb "red"     lw   2   t  'v_x'      ,\
     FILE_A    u  2:4  w   lines    lt  1   lc rgb "red"     lw   2   t  '-v_x'     ,\
     FILE_SG   u  1:2  w   points   pt  5   lc rgb "yellow"  ps   2   t  'Sol-Gio'

