from pylab import *

rc( "font" , size=20 )

# load single thread performance as reference

sizes = [ 512 , 1024 ]

threads = [ 2 , 4 , 8 , 16 ]

for N in sizes:

    data = loadtxt("../perf_results/perf_marvin_omp/perf_marvin_N%d_serial.dat" % N)
    runtime_serial = data[1]

    rt_omp = zeros( len(threads) )
    rt_df =  zeros( len(threads) )
    rt_ldf =  zeros( len(threads) )
    rt_ldfgb =  zeros( len(threads) )
    
    for i in xrange(len(threads)):
        data = loadtxt("../perf_results/perf_marvin_omp/perf_marvin_N%d_%d.dat" % (N,threads[i]))
        rt_omp[i] = min(data[:,1])
        data = loadtxt("../perf_results/perf_marvin_dataflow/perf_marvin_N%d_%d.dat" % (N,threads[i]))
        rt_df[i] = min(data[:,1])
        data = loadtxt("../perf_results/perf_marvin_local_dataflow/perf_marvin_N%d_%d.dat" % (N,threads[i]))
        rt_ldf[i] = min(data[:,1])
        data = loadtxt("../perf_results/perf_marvin_local_dataflow_gb/perf_marvin_N%d_%d.dat" % (N,threads[i]))
        rt_ldfgb[i] = min(data[:,1])

    figure()
    plot( threads , runtime_serial/rt_omp , 'o-' , label="OMP" )
    plot( threads , runtime_serial/rt_df , 'o-' , label="HPX df" )
    plot( threads , runtime_serial/rt_ldf , 'o-' , label="HPX ldf" )
    plot( threads , runtime_serial/rt_ldfgb , 'o-k' , label="HPX ldf gb" )

    title("System Size %dx%d" % (N,N))
    xlabel( "Threads" )
    ylabel( "Speedup" )
    legend( loc="upper left" )

    axis([1.5,16.5,1,14])

    figure()
    plot( threads , rt_omp , 'o-' , label="OMP" )
    plot( threads , rt_df , 'o-' , label="HPX df" )
    plot( threads , rt_ldf , 'o-' , label="HPX ldf" )
    plot( threads , rt_ldfgb , 'o-k' , label="HPX ldf gb" )

    title("System Size %dx%d" % (N,N))
    xlabel( "Threads" )
    ylabel( "Run Time (s)" )
    legend( loc="upper right" )

show()
