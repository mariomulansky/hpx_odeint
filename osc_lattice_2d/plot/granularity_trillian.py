from pylab import *

rc( "font" , size=20 )

N = 2048
t = 16

data_omp = loadtxt("../perf_results/perf_trillian_omp/perf_trillian_N%d_%d.dat" % (N,t))
data_df = loadtxt("../perf_results/perf_trillian_dataflow/perf_trillian_N%d_%d.dat" % (N,t))
data_ldf = loadtxt("../perf_results/perf_trillian_local_dataflow/perf_trillian_N%d_%d.dat" % (N,t))

plot( data_omp[:,0] , data_omp[:,1] , 'o-' , label="OMP" )
plot( data_df[1:,0] , data_df[1:,1] , 'o-' , label="HPX df" )
plot( data_ldf[:,0] , data_ldf[:,1] , 'o-' , label="HPX ldf" )

plot( data_ldf[:,0] , ones(len(data_ldf[:,0]))*min(data_ldf[:,1]) , '--k' )

title( "System Size %d, %d Threads" % (N,t) )

xlabel( "Granularity" )
ylabel( "Runtime (s)" )

legend( loc="upper right" )

show()
