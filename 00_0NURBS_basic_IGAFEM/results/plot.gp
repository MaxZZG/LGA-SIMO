      
      
      set term epslatex size 8cm,6cm color colortext lw 3			 	
      set out 'bar1d-p2e4.tex'

      set datafile separator ','
      set format x '$%3.1f$'                    
      set format y '$%3.2f$'                     
      set xtic auto                          
      set ytic 0.02                          
      set xlabel '$x$'
      set ylabel 'Displacement $u$'  rotate by 90

      set border linewidth 0.7
      set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 10
      set style line 2 linecolor rgb '#dd181f' linetype 2 linewidth 5

      plot    "bar1d1.csv" using 1:2 title 'Exact' with lines lt 1,\
              "bar1d2.csv" using 1:2 title 'IGA,p=2' with points
