from pylab import *

# load single thread performance as reference


rc( "font" , size=20 )
sizes = [ 1024 , 2048 ]

threads = [ 4 , 8 , 16 , 32 , 64]

for N in sizes:

    data = loadtxt("../perf_results/perf_trillian_omp/perf_trillian_N%d_serial.dat" % N)
    runtime_serial = data[1]

    rt_omp = zeros( len(threads) )
    rt_df =  zeros( len(threads) )
    rt_ldf =  zeros( len(threads) )
    
    for i in xrange(len(threads)):
        data = loadtxt("../perf_results/perf_trillian_omp/perf_trillian_N%d_%d.dat" % (N,threads[i]))
        rt_omp[i] = min(data[:,1])
        data = loadtxt("../perf_results/perf_trillian_dataflow/perf_trillian_N%d_%d.dat" % (N,threads[i]))
        rt_df[i] = min(data[:,1])
        data = loadtxt("../perf_results/perf_trillian_local_dataflow/perf_trillian_N%d_%d.dat" % (N,threads[i]))
        rt_ldf[i] = min(data[:,1])

    figure()
    plot( threads , runtime_serial/rt_omp , 'o-' , label="OMP" )
    plot( threads , runtime_serial/rt_df , 'o-' , label="HPX df" )
    plot( threads , runtime_serial/rt_ldf , 'o-' , label="HPX ldf" )

    title("System Size %dx%d" % (N,N))
    xlabel( "Threads" )
    ylabel( "Speedup" )
    legend( loc="upper left" )

show()
