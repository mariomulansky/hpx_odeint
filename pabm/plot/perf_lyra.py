from pylab import *

rc( "font" , size=20 )

cores = array([ 1 , 2 , 3 , 4 ])
rt3 = array([ 2.7588 , 1.401208 , 1.453303 , 0.866484 ])
rt4 = array([ 2.7588 , 1.401208 , 1.453303 , 0.866484 ])

plot( cores , rt[0]/rt , 'o-k' )

xlabel("Cores")
ylabel("Speed Up")

xticks( [1,2,3,4] , [1,2,3,4] )

axis( [0.7 , 4.5 , 0.5 , 4.1 ])

show()
