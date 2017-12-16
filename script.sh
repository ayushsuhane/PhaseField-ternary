arg1=$1
gnuplot -presist <<PLOT
set pm3d
unset surface
set view 0,0
set view equal xy
set format cb "%3.1f"
set palette defined (1.00 "red", 2.00 "green", 3.00 "blue", 4.00 "black")
#set cbtics ("Pearlite" 2.00 , "Ferrite" 2.50, "Austenite" 3.00 )
splot '$arg1' u 1:2:3
PLOT
