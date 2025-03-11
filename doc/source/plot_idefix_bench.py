import matplotlib.pyplot as plt
import json

def do_plot(title, bench_file, gpumodels):
    with open(bench_file, 'r') as f:
        benches = json.load(f)

    plt.figure()
    xmax=0
    ymax=0
    for gpumodel in gpumodels:
        select = [bench for bench in benches if bench['gpumodel'] == gpumodel][-1]

        xs = [r['nbgpu'] for r in select['results']]
        ys = [r['cell_updates'] for r in select['results']]
        plt.plot(xs, ys,'o-',label=gpumodel)
        xmax=max(xmax,max(xs))
        ymax=max(ymax,max(ys))

    plt.xscale("log", base=2)
    plt.ylim(0,ymax*1.1)
    plt.xlim(1,xmax*1.1)
    plt.legend()
    plt.xlabel("Number of GPUs/GCDs")
    plt.ylabel("Performance (cells / second / GPU)")
    plt.title(title)
